function sc_b_dynamics_xdot(t,x,m,params)


B_eci = interp1(params.t_vec,params.B_save,t)


max_moments = params.max_moments


p = x[1:3]
w = x[4:6]
# m = copy(u)

# p[1] = p[1]
# p[2] = p[2]
# p[3] = p[3]

#TODO: this isn't really the attitude, we just zero it out from ideal
# n_Q_b = dcm_from_p(p)
# b_Q_n = n_Q_b'
# B_body = b_Q_n*B_eci

B_body = dcm_from_p(p)'*B_eci

# w[1] = w[1]
# w[2] = w[2]
# w[3] = w[3]

# J = params.J
# inB_eciJ = params.inB_eciJ


# pdot = pdot_from_w(p,w)

# actual_m = max_moments .* tanh.(m)

# actual_tau= cross(max_moments .* tanh.(m),B_body)
# wdot = params.invJ*(actual_tau-cross(w,params.J*w))
#
# xdot = [pdot;wdot]
#
# dpdot_dp = [ ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)        (w[3]/2 - (p[2]*w[1])/2 + (p[1]*w[2])/2)      ((p[1]*w[3])/2 - (p[3]*w[1])/2 - w[2]/2) ;
#       ((p[2]*w[1])/2 - w[3]/2 - (p[1]*w[2])/2) ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)      (w[1]/2 - (p[3]*w[2])/2 + (p[2]*w[3])/2) ;
#       (w[2]/2 + (p[3]*w[1])/2 - (p[1]*w[3])/2)      ((p[3]*w[2])/2 - w[1]/2 - (p[2]*w[3])/2)   ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)]
#
# dpdot_dw = [ (p[1]^2/4 - p[2]^2/4 - p[3]^2/4 + 1/4)                 ((p[1]*p[2])/2 - p[3]/2)                 (p[2]/2 + (p[1]*p[3])/2);
#                (p[3]/2 + (p[1]*p[2])/2) (- p[1]^2/4 + p[2]^2/4 - p[3]^2/4 + 1/4)                 ((p[2]*p[3])/2 - p[1]/2);
#                ((p[1]*p[3])/2 - p[2]/2)                 (p[1]/2 + (p[2]*p[3])/2)  ( - p[1]^2/4 - p[2]^2/4 + p[3]^2/4 + 1/4)]
#
#
# # B_eci = B_eci;
# part1 = (-8*hat(p)*hat(B_eci) - 8*hat(hat(p)*B_eci) + 4*(1-p'*p)*hat(B_eci) + 8*hat(p)*B_eci*p')/(1+p'*p)^2
# part2 = -4*(8*hat(p)*hat(p)*B_eci - 4*(1-p'*p)*hat(p)*B_eci)*p'/((1+p'*p)^3)
# dwdot_dp = params.invJ*hat(max_moments .* tanh.(m))*(part1+part2)
#
# dwdot_dw = params.invJ*(hat(params.J*w) - hat(w)*params.J)
#
# A = [dpdot_dp dpdot_dw;
#      dwdot_dp dwdot_dw];
#
#
#
# taux = m[1]
# tauy = m[2]
# tauz = m[3]
#
# tauB = diagm(max_moments)*[ -1*(tanh(1*m[1])^2 - 1)      0   0;
#                            0       -1*(tanh(1*m[2])^2 - 1)   0;
#                            0        0   -1*(tanh(1*m[3])^2 - 1)]
#
#
# B = [zeros(3,3);
#      invJ*hat(-B_body)*tauB]
#


      return [pdot_from_w(p,w);params.invJ*(cross(params.max_moments .* tanh.(m),B_body)-cross(w,params.J*w))]
end


function sc_b_dynamics_full(t,x,m,params)


B_eci = interp1(params.t_vec,params.B_save,t)


max_moments = params.max_moments


p = x[1:3]
w = x[4:6]
# m = copy(u)

# p[1] = p[1]
# p[2] = p[2]
# p[3] = p[3]

#TODO: this isn't really the attitude, we just zero it out from ideal
# n_Q_b = dcm_from_p(p)
# b_Q_n = n_Q_b'
# B_body = b_Q_n*B_eci

B_body = dcm_from_p(p)'*B_eci

# w[1] = w[1]
# w[2] = w[2]
# w[3] = w[3]

# J = params.J
# inB_eciJ = params.inB_eciJ


pdot = pdot_from_w(p,w)

# actual_m = max_moments .* tanh.(m)

actual_tau= cross(max_moments .* tanh.(m),B_body)
wdot = params.invJ*(actual_tau-cross(w,params.J*w))

xdot = [pdot;wdot]

dpdot_dp = [ ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)        (w[3]/2 - (p[2]*w[1])/2 + (p[1]*w[2])/2)      ((p[1]*w[3])/2 - (p[3]*w[1])/2 - w[2]/2) ;
      ((p[2]*w[1])/2 - w[3]/2 - (p[1]*w[2])/2) ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)      (w[1]/2 - (p[3]*w[2])/2 + (p[2]*w[3])/2) ;
      (w[2]/2 + (p[3]*w[1])/2 - (p[1]*w[3])/2)      ((p[3]*w[2])/2 - w[1]/2 - (p[2]*w[3])/2)   ((p[1]*w[1])/2 + (p[2]*w[2])/2 + (p[3]*w[3])/2)]

dpdot_dw = [ (p[1]^2/4 - p[2]^2/4 - p[3]^2/4 + 1/4)                 ((p[1]*p[2])/2 - p[3]/2)                 (p[2]/2 + (p[1]*p[3])/2);
               (p[3]/2 + (p[1]*p[2])/2) (- p[1]^2/4 + p[2]^2/4 - p[3]^2/4 + 1/4)                 ((p[2]*p[3])/2 - p[1]/2);
               ((p[1]*p[3])/2 - p[2]/2)                 (p[1]/2 + (p[2]*p[3])/2)  ( - p[1]^2/4 - p[2]^2/4 + p[3]^2/4 + 1/4)]


# B_eci = B_eci;
part1 = (-8*hat(p)*hat(B_eci) - 8*hat(hat(p)*B_eci) + 4*(1-p'*p)*hat(B_eci) + 8*hat(p)*B_eci*p')/(1+p'*p)^2
part2 = -4*(8*hat(p)*hat(p)*B_eci - 4*(1-p'*p)*hat(p)*B_eci)*p'/((1+p'*p)^3)
dwdot_dp = params.invJ*hat(max_moments .* tanh.(m))*(part1+part2)

dwdot_dw = params.invJ*(hat(params.J*w) - hat(w)*params.J)

A = [dpdot_dp dpdot_dw;
     dwdot_dp dwdot_dw];



# taux = m[1]
# tauy = m[2]
# tauz = m[3]

tauB = diagm(max_moments)*[ -1*(tanh(1*m[1])^2 - 1)      0   0;
                           0       -1*(tanh(1*m[2])^2 - 1)   0;
                           0        0   -1*(tanh(1*m[3])^2 - 1)]


B = [zeros(3,3);
     invJ*hat(-B_body)*tauB]



      return A, B
end
## here we forward diff the dynamics
t = 403.4
u = randn(3)

x_test = randn(6)


fwd_f(x) = sc_b_dynamics_xdot(t,x,u,params)
J = x -> ForwardDiff.jacobian(fwd_f,x)


xdot,A,B = sc_b_dynamics_full(t,x_test,u,params)

Afd = J(x_test)


fwd_u(u) = sc_b_dynamics_xdot(t,x_test,u,params)

J_2 = u -> ForwardDiff.jacobian(fwd_u,u)

Bfd = J_2(u)
