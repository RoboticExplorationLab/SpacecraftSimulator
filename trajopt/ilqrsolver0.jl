function iLQRsimple_B(x0, xg, utraj, Q, R, Qf, dt,params)
# %#codegen
# %iLQR Trajectory Optimization
dJtol = params.dJ_tol

# % unpack stuff
# eps = params.eps;
# Q_scale = params.Q_scale;
# end_bump = params.end_bump;
Nx = length(x0)
Nu = size(utraj,1)
N = size(utraj,2)+1

xtraj = zeros(Nx, N)
xtraj = convert(Array{Float32,2},xtraj)
xtraj[:,1] = x0

A = convert(Array{Float32,3},zeros(Nx, Nx, N-1))
B = convert(Array{Float32,3},zeros(Nx, Nu, N-1))

# %First simulate with utraj0 to get initial matrices
J = 0.0
for k = 1:(N-1)
    J = J + (1/2)*(xtraj[:,k]-xg)'*Q*(xtraj[:,k]-xg) + (1/2)*utraj[:,k]'*R*utraj[:,k]
    xtraj[:,k+1], A[:,:,k], B[:,:,k] = rk4step(xtraj[:,k],utraj[:,k], dt,(k-1)*dt,params)
end
J = J+ (xtraj[:,N]-xg)'*Qf*(xtraj[:,N]-xg)

# @infiltrate
# error()

# % Set up backwards pass matrices
K = convert(Array{Float32,3},zeros(Nu,Nx,N-1))
l = convert(Array{Float32,2},zeros(Nu,N-1))


# iter = 0
max_iters = 50
for iter = 1:max_iters

    # iter = iter + 1;

    # %Set up backwards LQR pass
    S = Qf;
    s = Qf*(xtraj[:,N]-xg);

    for k = (N-1):-1:1

        # %Calculate cost gradients for this time step
        q = Q*(xtraj[:,k]-xg)
        r = R*utraj[:,k]

        # % inverses come in here
        ilqrinv = inv(convert(Array{Float32,2},(R + B[:,:,k]'*S*B[:,:,k])))
        # l[:,k] = (R + B[:,:,k]'*S*B[:,:,k])\(r + B[:,:,k]'*s)
        # K[:,:,k] = (R + B[:,:,k]'*S*B[:,:,k])\(B[:,:,k]'*S*A[:,:,k])
        l[:,k] = ilqrinv*(r + B[:,:,k]'*s)
        K[:,:,k] = ilqrinv*(B[:,:,k]'*S*A[:,:,k])

        # @infiltrate
        # error()
        # %Calculate new S and s
        Snew = Q + K[:,:,k]'*R*K[:,:,k] + (A[:,:,k]-B[:,:,k]*K[:,:,k])'*S*(A[:,:,k]-B[:,:,k]*K[:,:,k])
        snew = q - K[:,:,k]'*r + K[:,:,k]'*R*l[:,k] + (A[:,:,k]-B[:,:,k]*K[:,:,k])'*(s - S*B[:,:,k]*l[:,k])
        S = Snew
        s = snew

    end

    # %Now do forward pass line search with new l and K
    unew = convert(Array{Float32,2},zeros(Nu,N-1))
    xnew = convert(Array{Float32,2},zeros(Nx,N))
    xnew[:,1] = xtraj[:,1];
    alpha = params.alpha0;
    Jnew = J+1;
    # while Jnew > J
    for line_i = 1:6
        Jnew = 0.0

        # @infiltrate
        for k = 1:N-1
            unew[:,k] = utraj[:,k] - alpha*l[:,k] - K[:,:,k]*(xnew[:,k]-xtraj[:,k])
            xnew[:,k+1], A[:,:,k], B[:,:,k] = rk4step(xnew[:,k],unew[:,k], dt,(k-1)*dt,params)

            Jnew = Jnew + (1/2)*(xnew[:,k]-xg)'*Q*(xnew[:,k]-xg) + (1/2)*unew[:,k]'*R*unew[:,k]

        end
        Jnew = Jnew + (1/2)*(xnew[:,N]-xg)'*Qf*(xnew[:,N]-xg)

        # @infiltrate
        # error()
        if Jnew<J # this also pulls the linesearch back if hasnan(xnew)
            break
        else
            alpha = (1/2)*alpha
        end

    end

    dJ = abs(J - Jnew)
    if dJ<dJtol
        break
    end


    if rem((iter-1),4)==0
        println("iter          alpha      maxL    Cost")
    end
    maxL = round(maximum(vec(l)),digits = 3)
    J = Jnew
    J_display = round(J,digits = 3)
    alpha_display = round(alpha,digits = 3)
    println("$iter          $alpha_display      $maxL    $J_display")


    xtraj = xnew
    utraj = unew
end

return xtraj, utraj, K
end


#     for iii = 1:N
#         angle_error_hist{iter}(iii) = norm(phi_from_p(xtraj(1:3,iii)));
#     end
#
# %     if abs(J-Jnew)<dJtol
# %         break
# %     else
#     J = Jnew
# %     end

# @show iter

# if rem((iter-1),4)==0
#     println("iter          alpha      maxL    Cost")
# end
# ll = max(max(abs(l)));
#
#
# disp([iter 2*alpha ll J])
# end
#
# end

# function [x1, A, B] = rkstep(x0,u0,dt,t,params)
#     %Explicit midpoint step from x_k to x_{k+1}
#     Nx = length(x0);
#     Nu = length(u0);
#
# %     [xdot1, dxdot1] = pendulum_dynamics(0,x0,u0);
# %     [xdot2, dxdot2] = pendulum_dynamics(0,x0+.5*dt*xdot1,u0);
#     [xdot1, A1, B1] = sc_b_dynamics(t,x0,u0,params);
#     [xdot2, A2, B2] = sc_b_dynamics(t+.5*dt,x0+.5*dt*xdot1,u0,params);
#     x1 = x0 + dt*xdot2;
#
# %     A1 = dxdot1(:,(1:Nx));
# %     B1 = dxdot1(:,Nx+(1:Nu));
# %
# %     A2 = dxdot2(:,(1:Nx));
# %     B2 = dxdot2(:,Nx+(1:Nu));
#
#     A = eye(Nx) + dt*A2 + 0.5*dt*dt*A2*A1;
#     B = dt*B2 + 0.5*dt*dt*A2*B1;
#
# %     disp('ye')
# end

function rk4step(x0,u0,dt,t,params)

    Nx = length(x0)
# %     Nu = length(u0);

    x1 = x0

    xdot1, A1, B1 = sc_b_dynamics(t,x0,u0,params)
    k1 = xdot1*dt

    x2 = x1 + .5*k1
    xdot2, A2, B2 = sc_b_dynamics(t+dt*.5,x2,u0,params)
    k2 = dt*xdot2;

    x3 = x1 + .5*k2
    xdot3, A3, B3 = sc_b_dynamics(t+dt*.5,x3,u0,params)
    k3 = dt*xdot3

    x4 = x1 + k3
    xdot4, A4, B4 = sc_b_dynamics(t+dt,x4,u0,params)
    k4 = dt*xdot4

    x_tp1 = x1 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    # midpoint method
    # A_d = eye(Nx) + dt*A2 + 0.5*dt*dt*A2*A1
    # B_d = dt*B2 + 0.5*dt*dt*A2*B1

    # RK4 method
    #A_d
    dk1_dx1 = dt*A1
    dx2_dx1 = eye(Nx) + .5*dk1_dx1
    dk2_dx1 = dt*A2*dx2_dx1
    dx3_dx1 = eye(Nx) + .5*dk2_dx1
    dk3_dx1 = dt*A3*dx3_dx1
    dx4_dx1 = eye(Nx) + dk3_dx1
    dk4_dx1 = dt*A4*dx4_dx1
    A_d = eye(Nx) + (1/6)*(dk1_dx1 + 2*dk2_dx1 + 2*dk3_dx1 + dk4_dx1)

    # B_d
    dk1_du = dt*B1
    dx2_du = .5*dk1_du
    dk2_du = dt*A2*dx2_du + dt*B2
    dx3_du = .5*dk2_du
    dk3_du = dt*A3*dx3_du + dt*B3
    dx4_du = dk3_du
    dk4_du = dt*A4*dx4_du + dt*B4
    B_d = (1/6)*(dk1_du + 2*dk2_du + 2*dk3_du + dk4_du)


    return x_tp1, A_d, B_d
end
