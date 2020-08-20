function quad(Q::Array{Float32,2},
              x::Union{Array{Float32,1},Array{Float32,2}})
    """returns the quadratic xᵀQx"""
    return x'*Q*x
end

function cost(Q::Array{Float32,2},R::Array{Float32,2},
              x::Array{Float32,1},u::Array{Float32,1})
    """returns quadratic only cost for given Q,R,x, and u"""
    return .5*x'*Q*x + .5*u'*R*u
end



function iLQRsimple_B3(x0::Array{Float32,1},
                       xg::Array{Float32,1},
                       Q::Array{Float32,2},
                       R::Array{Float32,2},
                       Qf::Array{Float32,2},
                       N::Int,
                       dt::Float32,
                       params::NamedTuple)::Tuple{Array{Array{Float32,1},1},
                                                  Array{Array{Float32,1},1},
                                                  Array{Array{Float32,2},1}}
"""ILQR with all vectors of vectors/mats, and Float32

Args:
    x0       : initial condition [ᴺpᴮ;ω] ∈ R⁶
    xg       : goal state
    Q        : state cost quadratic term
    R        : control cost quadratic term
    Qf       : final state cost quadratic term
    N        : number of knot points
    dt       : timestep between knot points
    params   : named tuple with spacecraft info

Returns:
    xtraj    : vector of vectors for state history
    utraj    : vector of vectors for control history
    K        : vector of gain matrices
"""

# state and control dimensions
Nx = 6
Nu = 3

# initial trajectory is initial conditions the whole time
xtraj = fill(x0,N)
utraj = fill(zeros(Float32,Nu),N-1)

# cost corresponding to no control
J = (N-1)*quad(Q,x0-xg) + quad(Qf,(xtraj[N]-xg))

# allocate K and l
K = fill(zeros(Float32,Nu,Nx),N-1)
l = fill(zeros(Float32,Nu),N-1)

# allocate the new states and controls
xnew = fill(zeros(Float32,Nx),N)
unew = fill(zeros(Float32,Nu),N-1)

# main loop
for iter = 1:50

    # cost to go matrices at the end of the trajectory
    S = Qf
    s = Qf*(xtraj[N]-xg)

    # backwards pass
    for k = (N-1):-1:1

        # calculate cost gradients for this time step
        q = Q*(xtraj[k]-xg)
        r = R*utraj[k]

        # jacobians
        Ak, Bk = rk4step_jacobians(xtraj[k],utraj[k], dt,(k-1)*dt,params)

        # linear solve
        LHinv = inv(R + Bk'*S*Bk)
        l[k] = LHinv*(r + Bk'*s)
        K[k] = LHinv*(Bk'*S*Ak)

        # update
        Snew = Q + K[k]'*R*K[k] + quad(S,(Ak-Bk*K[k]))
        snew = q - K[k]'*r + K[k]'*R*l[k] + (Ak-Bk*K[k])'*(s - S*Bk*l[k])

        # update S's
        S = copy(Snew)
        s = copy(snew)

    end

    # initial conditions
    xnew[1] = x0

    # learning rate
    alpha = 1.0

    # line search
    Jnew = 0.0
    for line_i = 1:6
        Jnew = 0.0

        # rollout the dynamics
        for k = 1:N-1
            unew[k] = utraj[k] - alpha*l[k] - K[k]*(xnew[k]-xtraj[k])
            xnew[k+1]  = rk4step_xdot(xnew[k],unew[k], dt,(k-1)*dt,params)
            Jnew += cost(Q,R,xnew[k],unew[k])
        end
        Jnew = Jnew + (1/2)*(xnew[N]-xg)'*Qf*(xnew[N]-xg)

        # if the new cost is lower, we keep it
        if Jnew<J
            break
        else# this also pulls the linesearch back if hasnan(xnew)
            alpha = (1/2)*alpha
        end

    end

    # update trajectory and control history
    xtraj = copy(xnew)
    utraj = copy(unew)

    # termination criteria
    if abs(J - Jnew)<params.dJ_tol
        break
    end


    # ----------------------------output stuff-----------------------------
    if rem((iter-1),4)==0
        println("iter          alpha      maxL    Cost")
    end
    maxL = round(maximum(vec(mat_from_vec(l))),digits = 3)
    J = Jnew
    J_display = round(J,digits = 3)
    alpha_display = round(alpha,digits = 3)
    println("$iter          $alpha_display      $maxL    $J_display")



end

return xtraj, utraj, K
end


function rk4step_xdot(x1,u0,dt,t,params)

    # Nx = length(x0)

    # x1 = x0

    xdot1 = sc_b_dynamics_xdot(t,x0,u0,params)
    k1 = xdot1*dt

    x2 = x1 + .5*k1
    xdot2= sc_b_dynamics_xdot(t+dt*.5,x2,u0,params)
    k2 = dt*xdot2;

    x3 = x1 + .5*k2
    xdot3 = sc_b_dynamics_xdot(t+dt*.5,x3,u0,params)
    k3 = dt*xdot3

    x4 = x1 + k3
    xdot4= sc_b_dynamics_xdot(t+dt,x4,u0,params)
    k4 = dt*xdot4

    x_tp1 = x1 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    x_tp1 = convert(Array{Float32,1},x_tp1)

    return x_tp1
end
#
function rk4step_jacobians(x0,u0,dt,t,params)

    Nx = length(x0)

    x1 = x0

    xdot1, A1, B1 = sc_b_dynamics(t,x0,u0,params)
    k1 = xdot1*dt

    x2 = x1 + .5*k1
    xdot2, A2, B2 = sc_b_dynamics(t+dt*.5,x2,u0,params)
    k2 = dt*xdot2;
    #
    # x3 = x1 + .5*k2
    # xdot3, A3, B3 = sc_b_dynamics(t+dt*.5,x3,u0,params)
    # k3 = dt*xdot3
    #
    # x4 = x1 + k3
    # xdot4, A4, B4 = sc_b_dynamics(t+dt,x4,u0,params)
    # k4 = dt*xdot4
    #
    # x_tp1 = x1 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    # midpoint method
    A_d = eye(Nx) + dt*A2 + 0.5*dt*dt*A2*A1
    B_d = dt*B2 + 0.5*dt*dt*A2*B1

    # RK4 method
    #A_d
    # dk1_dx1 = dt*A1
    # dx2_dx1 = eye(Nx) + .5*dk1_dx1
    # dk2_dx1 = dt*A2*dx2_dx1
    # dx3_dx1 = eye(Nx) + .5*dk2_dx1
    # dk3_dx1 = dt*A3*dx3_dx1
    # dx4_dx1 = eye(Nx) + dk3_dx1
    # dk4_dx1 = dt*A4*dx4_dx1
    # A_d = eye(Nx) + (1/6)*(dk1_dx1 + 2*dk2_dx1 + 2*dk3_dx1 + dk4_dx1)
    #
    # # B_d
    # dk1_du = dt*B1
    # dx2_du = .5*dk1_du
    # dk2_du = dt*A2*dx2_du + dt*B2
    # dx3_du = .5*dk2_du
    # dk3_du = dt*A3*dx3_du + dt*B3
    # dx4_du = dk3_du
    # dk4_du = dt*A4*dx4_du + dt*B4
    # B_d = (1/6)*(dk1_du + 2*dk2_du + 2*dk3_du + dk4_du)


    A_d   = convert(Array{Float32,2},A_d)
    B_d   = convert(Array{Float32,2},B_d)
    return A_d, B_d
end












# scope test
function testit()


    a = 3.0

    for i = 1:3
        a += 1.0
    end

    @show a

end
