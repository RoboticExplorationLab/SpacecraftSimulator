function sc_b_dynamics(t,x,u,params)

# @infiltrate
# error()
# B_eci = interp1(params.t_vec,params.B_save,t)

T = 50*60
B_eci = 3e-5*[sin(2*pi*t/T + pi/2);sin(2*pi*t/T - pi/3);sin(2*pi*t/T + 3*pi/2)]

max_moments = params.max_moments
# x = x[:]
p = x[1:3]
w = x[4:6]
m = copy(u)

p1 = p[1]
p2 = p[2]
p3 = p[3]

#TODO: this isn't really the attitude, we just zero it out from ideal
n_Q_b = dcm_from_p(p)
b_Q_n = n_Q_b'
B_body = b_Q_n*B_eci

wx = w[1]
wy = w[2]
wz = w[3]

J = params.J
invJ = params.invJ


pdot = pdot_from_w(p,w)

actual_m = max_moments .* tanh.(m)

actual_tau= cross(actual_m,B_body)
wdot = invJ*(actual_tau-cross(w,J*w))

xdot = [pdot;wdot]

dpdot_dp = [ ((p1*wx)/2 + (p2*wy)/2 + (p3*wz)/2)        (wz/2 - (p2*wx)/2 + (p1*wy)/2)      ((p1*wz)/2 - (p3*wx)/2 - wy/2) ;
      ((p2*wx)/2 - wz/2 - (p1*wy)/2) ((p1*wx)/2 + (p2*wy)/2 + (p3*wz)/2)      (wx/2 - (p3*wy)/2 + (p2*wz)/2) ;
      (wy/2 + (p3*wx)/2 - (p1*wz)/2)      ((p3*wy)/2 - wx/2 - (p2*wz)/2)   ((p1*wx)/2 + (p2*wy)/2 + (p3*wz)/2)]

dpdot_dw = [ (p1^2/4 - p2^2/4 - p3^2/4 + 1/4)                 ((p1*p2)/2 - p3/2)                 (p2/2 + (p1*p3)/2);
               (p3/2 + (p1*p2)/2) (- p1^2/4 + p2^2/4 - p3^2/4 + 1/4)                 ((p2*p3)/2 - p1/2);
               ((p1*p3)/2 - p2/2)                 (p1/2 + (p2*p3)/2)  ( - p1^2/4 - p2^2/4 + p3^2/4 + 1/4)]


v = B_eci;
part1 = (-8*hat(p)*hat(v) - 8*hat(hat(p)*v) + 4*(1-p'*p)*hat(v) + 8*hat(p)*v*p')/(1+p'*p)^2
part2 = -4*(8*hat(p)*hat(p)*v - 4*(1-p'*p)*hat(p)*v)*p'/((1+p'*p)^3)
dwdot_dp = invJ*hat(actual_m)*(part1+part2)

dwdot_dw = invJ*(hat(J*w) - hat(w)*J)

A = [dpdot_dp dpdot_dw;
     dwdot_dp dwdot_dw];



taux = m[1]
tauy = m[2]
tauz = m[3]

tauB = diagm(max_moments)*[ -1*(tanh(1*taux)^2 - 1)      0   0;
                           0       -1*(tanh(1*tauy)^2 - 1)   0;
                           0        0   -1*(tanh(1*tauz)^2 - 1)]


B = [zeros(3,3);
     invJ*hat(-B_body)*tauB]

# % dxdot = [A,B];

      return xdot, A, B
end


function sc_b_dynamics_xdot(t,x,u,params)

      # u = interp1(params.u_t,uhist,t)
      xdot, A, B = sc_b_dynamics(t,x,u,params)

      return xdot
end
# function sc_b_dynamics_jacobians(t,x,u,params)
#
#       # u = interp1(params.u_t,uhist,t)
#       xdot, A, B = sc_b_dynamics(t,x,u,params)
#
#       return A, B
# end

function rk4(f,t_n,x_n,u,params,dt)
    # standard rk4

    k1 = dt*f(t_n,x_n,u,params)
    k2 = dt*f(t_n+dt/2, x_n+k1/2,u,params)
    k3 = dt*f(t_n+dt/2, x_n+k2/2,u,params)
    k4 = dt*f(t_n+dt, x_n+k3,u,params)

    x_np1 = (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

        return x_np1
end


## here we forward diff the dynamics
# t = 403.4
# u = randn(3)
#
# x_test = randn(6)
#
#
# fwd_f(x) = sc_b_dynamics2(t,x,u,params)
# J = x -> ForwardDiff.jacobian(fwd_f,x)
#
#
# xdot,A,B = sc_b_dynamics(t,x_test,u,params)
#
# Afd = J(x_test)
#
#
# fwd_u(u) = sc_b_dynamics2(t,x_test,u,params)
#
# J_2 = u -> ForwardDiff.jacobian(fwd_u,u)
#
# Bfd = J_2(u)
