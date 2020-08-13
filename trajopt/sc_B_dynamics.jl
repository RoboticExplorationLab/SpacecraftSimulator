function sc_b_dynamics(t,x,u,params)

# % B_eci_fx = @(t) .001*expm(hat([0;1/sqrt(2);1/sqrt(2)]*t*.1*2*pi/500))*[1;0;0];

# % B_eci = B_eci_fx(t);


# B_1 = interp1(params.t_vec,params.B(:,1),t);
# B_2 = interp1(params.t_vec,params.B(:,2),t);
# B_3 = interp1(params.t_vec,params.B(:,3),t);

# B_eci = [B_1;B_2;B_3];
B_eci = interp1(params.t_vec,params.B_save,t)

# u_scale = params.u_scale;
#
# u = u/u_scale;

# 1 = params.1;
# tau_mag = params.tau_mag;



max_moments = params.max_moments
# x = x[:]
p = x[1:3]
w = x[4:6]
m = copy(u)

p1 = p[1]
p2 = p[2]
p3 = p[3]

n_Q_b = dcm_from_p(p)
b_Q_n = n_Q_b'
B_body = b_Q_n*B_eci

wx = w[1]
wy = w[2]
wz = w[3]

J = params.J
invJ = params.invJ
# J = [1.959e-4 2016.333e-9 269.176e-9;...
#            2016.333e-9 1.999e-4 2318.659e-9;...
#            269.176e-9 2318.659e-9 1.064e-4];
# invJ = 1e3*[5.105190043647774  -0.051357738055630  -0.011796180015288;...
#   -0.051357738055630   5.004282688631113  -0.108922940075563;...
#   -0.011796180015288  -0.108922940075563   9.400899721840833];

pdot = pdot_from_w(p,w)

actual_m = max_moments .* tanh.(m)

# % if max(actual_m)> tau_mag
# %     disp('fuck')
# % end

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

tauB = diagm(max_moments)*[ -1*(tanh(1*taux)^2 - 1)                   0                    0;
                                         0   -1*(tanh(1*tauy)^2 - 1)               0;
                                         0                                         0   -1*(tanh(1*tauz)^2 - 1)]


B = [zeros(3,3);
     invJ*hat(-B_body)*tauB]

# % dxdot = [A,B];

      return xdot, A, B
end
