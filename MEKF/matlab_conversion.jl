using LinearAlgebra

mutable struct EKF_struct
    mu::Array{Array{Float64,1},1}
    Sigma::Array{Array{Float64,2},1}
end

function dynODE(t_n,x_n,tau,params)

    J = params.J
    invJ = params.invJ

    # x_n = x_n[:];
    q = x_n[1:4]
    w = x_n[5:7]
    # beta = x_n(4:6);

    q_dot = .5*qdot(q,[w;0])

    wdot = invJ*(tau - cross(w,J*w))

    beta_dot = zeros(3)

    Xdot = [q_dot;
            wdot;
            beta_dot]

    return Xdot

end

function g(x,r1n, r2n)

    x = x[:];
    q = x[1:4];
    w = x[5:7];

    beta = x[8:10];

    gyro = w+beta;

    n_Q_b = dcm_from_q(q);
    b_Q_n = n_Q_b';

    r1b = b_Q_n*r1n;
    r2b = b_Q_n*r2n;

    y = [gyro;r1b;r2b];

    return y

end

function add_noise_to_y(y,sqrtR)
# add measurement noise t
# y = y[:];
gyro = y[1:3];
r1b = y[4:6];
r2b = y[7:9];

noisy_gyro = gyro + fast_mvnrnd(sqrtR[1:3,1:3]);
noisy_r1b = S03_noise(r1b,sqrtR[4:6,4:6]);
noisy_r2b = S03_noise(r2b,sqrtR[7:9,7:9]);

noisy_y = [noisy_gyro;noisy_r1b;noisy_r2b];

return noisy_y

end


function S03_noise(vector, sqrt_covariance)
    # here we add noise S03 style

    # noise axis angle vector
    phi_noise = fast_mvnrnd(sqrt_covariance);


    # apply the noise rotation to the vector
    noisy_vector = skew_expm(hat(phi_noise))*vector;

    return noisy_vector
end



function ekf_predict(mu,Sigma,tau,params,dt,Q)

    # tn = 0;

    # unpack state
    # mu = mu[:];
    q = mu[1:4];
    w = mu[5:7];
    beta = mu[8:10];

    # inertia
    J = params.J;
    invJ = params.invJ;

    # predict using RK4
    # mu_predict = rk4(dynODE,tn,mu,tau,params,dt);
    q_next = qdot(q,q_from_phi(w*dt));
    wdot = invJ*(tau-cross(w,J*w));
    w_next = w + wdot*dt;
    beta_next = beta;

    mu_predict = [q_next;w_next;beta_next];

    # jacobians
    dphi_phi = I- hat(w*dt);
    dphi_dw = .5*dt*eye(3);
    dw_w = (I + dt*invJ*(hat(J*w) - hat(w)*J));
    dbeta_beta = eye(3);

    # dynamics jacobian (d_f/d_state)
    A = [dphi_phi        dphi_dw          zeros(3,3);
         zeros(3,3)      dw_w                zeros(3,3);
         zeros(3,3)      zeros(3,3)          dbeta_beta];

    # covariance prediction
    Sigma_predict = A*Sigma*A' + Q;

    return mu_predict, Sigma_predict
end


function ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R)

    # predicted quaternion
    q_predict = mu_predict[1:4]

    # dcm
    n_Q_bt = dcm_from_q(q_predict)

    # use the predicted attitude to generate predicted measurements
    r1b_t = (n_Q_bt)'*r1n
    r2b_t = (n_Q_bt)'*r2n


    # predicted measurment
    yhat = g(mu_predict,r1n, r2n)

    # measurement difference
    z = yt - yhat

    # measurement Jacobian
    C = [zeros(3,3)     I   I  ;
         2*hat(r1b_t) zeros(3,3) zeros(3,3);
         2*hat(r2b_t) zeros(3,3) zeros(3,3)]

    # covariance of the innovation
    S = C*Sigma_predict*C' + R

    # kalman gain
    # L = Sigma_predict*C'*inv(S)
    L = Sigma_predict*C'/S

    # innovation
    delta_mu = L*z

    # update quaternion multiplicatively
    q_update = qdot(q_predict,q_from_g(delta_mu[1:3]))

    # update angular velocity and bias additively
    w_update = mu_predict[5:7] + delta_mu[4:6];
    beta_update = mu_predict[8:10] + delta_mu[7:9];

    # updated mu
    mu_update = [ q_update ;
                  w_update;
                  beta_update]

    # Sigma update
    Sigma_update = (eye(9) - L*C)*Sigma_predict;

    return mu_update, Sigma_update
end

function MEKF(mu,Sigma,tau,yt,Q,R,dt,r1n,r2n,params)

# predict
mu_predict,Sigma_predict = ekf_predict(mu,Sigma,tau,params,dt,Q);

# innovate/update
mu_update,Sigma_update = ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R);

    return mu_update, Sigma_update
end


function triad_equal_error(r1,r2,b1,b2)
    # triad but with the errors spread equally between the two vectors
    # r1 and r2 are expressed in the inertial frame
    # b1 and b2 are the same vectors expressed in the body frame

    # these should be proportional to the inverse square of the covariance
    # for each measurement
    a1 = .5
    a2 = .5

    r3 = normalize(cross(r1,r2))
    b3 = normalize(cross(b1,b2))

    # this is from the FOAM method
    lambda_max = sqrt(a1^2 + a2^2 + 2*a1*a2*( dot(b1,b2)*dot(r1,r2)
                 + norm(cross(b1,b2))*norm(cross(r1,r2))  ))


    n_Q_b = (b3*r3' + (a1/lambda_max)*(b1*r1' + cross(b1,b3)*cross(r1,r3)')
             + (a2/lambda_max)*(b2*r2' + cross(b2,b3)*cross(r2,r3)'))'

    # get quaternion
    n_q_b = q_from_dcm(n_Q_b);

    return n_q_b, n_Q_b
end


function q_angle_error(n_q_b,n_q_b_est)
    # gets the error between the estimated and true attitude

    # error quaternion b_q_b_est
    errorq = q_shorter(qdot(qconj(n_q_b),n_q_b_est));


    # convert to phi
    errorphi = phi_from_q(errorq);

    # angle of the axis angle vector is the norm
    error_angle = norm(errorphi);

    return error_angle
end

function rk4(ODE,tn,xn,u,sc,h)
    # rk4 for a single step

    xn = xn[:]

    k1 = h*ODE(tn,xn,u,sc);
    k2 = h*ODE(tn + h/2,xn + k1/2,u,sc);
    k3 = h*ODE(tn + h/2,xn + k2/2,u,sc);
    k4 = h*ODE(tn + h,xn + k3,u,sc);

    y_np1 = xn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    return y_np1


end


# %% Sim
function run_MEKF()
# % inertia
J = diagm([1.0;2;3]);
invJ = inv(J);
# params.J = J; params.invJ = invJ;
params = (J=J,invJ = invJ)

# % initial conditions
q0 = randq()
w0 = deg2rad(5)*normalize(randn(3));
beta0 = deg2rad.([1;-1;.5]);

# % two known vectors expressed in ECI
r1n = normalize(randn(3));
r2n = normalize(randn(3));

# % time stuff
dt = .05;
tf = 2000.0;
t_vec = 0:dt:tf;

# % pre-allocate
state_hist = zeros(10,length(t_vec));
meas_hist = zeros(9,length(t_vec));
state_hist[1:4,1] = q0;
state_hist[5:7,1] = w0;
state_hist[8:10,1] = beta0;
meas_hist[:,1] = g([q0;w0;beta0],r1n,r2n);
q_est_hist = zeros(4,length(t_vec));
q_est_hist[:,1] = q0;
angle_error = zeros(length(t_vec));
ekf_angle_error = zeros(length(t_vec));

# % Process noise covariance
Q = zeros(10,10);
Q[8:10,8:10] = .0000001*eye(3);
Q_ekf = .0000000001*eye(9);
Q_ekf[7:9,7:9] = Q[8:10,8:10];

sqrtQ = sqrt(Q)
# % Sensor noise covariance
R_gyro = (deg2rad(1))^2*eye(3);
R_sun_sensor = (deg2rad(5))^2*eye(3);
R_magnetometer = (deg2rad(5))^2*eye(3);
R = zeros(9,9);
R[1:3,1:3] = R_gyro;
R[4:6,4:6] = R_sun_sensor;
R[7:9,7:9] = R_magnetometer;

sqrtR = sqrt(R)

# % external torque on spacecraft
tau = zeros(3);

mekf_time = 0.0

# % MEKF initialize
EKF = EKF_struct(copy(fill(zeros(10),length(t_vec))),copy(fill(zeros(9,9),length(t_vec))))

# mu = state_hist;
# Sigma = cell(length(t_vec),1);
EKF.mu[1] = state_hist[:,1]
EKF.Sigma[1] = .01*eye(9);

# % main loop
@showprogress "Simulating..." for kk = 1:(length(t_vec)-1)

    # % true dynamics
    # [t,Y] = ode45(@(t_n,x_n) dynODE(t_n,x_n,tau,params), [0,dt],state_hist(:,kk));
    # state_hist[:,kk+1] = rk4(dynODE,t_vec[kk],state_hist[:,kk],tau,params,dt) + mvnrnd(zeros(10),Q)
    state_hist[:,kk+1] = rk4(dynODE,t_vec[kk],state_hist[:,kk],tau,params,dt) + fast_mvnrnd(sqrtQ)
    # state_hist(:,kk+1) = Y(end,:)' + mvnrnd(zeros(10,1),Q)';
    state_hist[1:4,kk+1] = normalize(state_hist[1:4,kk+1]);

    # % true measurement
    meas_hist[:,kk+1] = add_noise_to_y(g(state_hist[:,kk+1],r1n,r2n),sqrtR);

    # % triad on the measurements
    r1b = meas_hist[4:6,kk+1];
    r2b = meas_hist[7:9,kk+1];
    n_q_b_est, n_Q_b_est = triad_equal_error(r1n,r2n,r1b,r2b);
    q_est_hist[:,kk+1] = n_q_b_est;

    # % triad angle error
    angle_error[kk+1] = q_angle_error(state_hist[1:4,kk+1],n_q_b_est);

# %     if rad2deg(angle_error(kk+1)) > 150
# %         n_Q_b = dcm_from_q(state_hist(1:4,kk+1));
# %         b_Q_n = n_Q_b';
# %         r1b = meas_hist(4:6,kk+1);
# %         r2b = meas_hist(7:9,kk+1);
# %
# %         disp('he')
# %     end
# %
    # MEKF
    t1 = time()
    EKF.mu[kk+1], EKF.Sigma[kk+1] = MEKF(EKF.mu[kk],EKF.Sigma[kk],tau,meas_hist[:,kk+1],Q_ekf,R,dt,r1n,r2n,params);
    mekf_time += time()-t1
    # % MEKF angle error
    ekf_angle_error[kk+1] = q_angle_error(state_hist[1:4,kk+1],EKF.mu[kk+1][1:4]);
end

# display(plot(t_vec,rad2deg.(ekf_angle_error),label = "MEKF"))
# display(plot!(t_vec,rad2deg.(angle_error),label = "Triad"))

println("Done!")
@show mekf_time

end


run_MEKF()
# %% plotting
# q_hist = state_hist(1:4,:);
# omega_hist = state_hist(5:7,:);
# beta_hist = state_hist(8:10,:);
#
# % figure
# % hold on
# % title('True \omega')
# % plot(t_vec,omega_hist(1,:))
# % plot(t_vec,omega_hist(2,:))
# % plot(t_vec,omega_hist(3,:))
# % hold off
#
# % figure
# % hold on
# % title('True quaternion')
# % plot(t_vec,q_hist(1,:))
# % plot(t_vec,q_hist(2,:))
# % plot(t_vec,q_hist(3,:))
# % plot(t_vec,q_hist(4,:))
# % hold off
#
# figure
# hold on
# title('est quaternion')
# plot(t_vec,q_est_hist(1,:))
# plot(t_vec,q_est_hist(2,:))
# plot(t_vec,q_est_hist(3,:))
# plot(t_vec,q_est_hist(4,:))
# hold off
#
#
# figure
# hold on
# title('Bias')
# plot(t_vec,rad2deg(beta_hist(1,:)))
# plot(t_vec,rad2deg(beta_hist(2,:)))
# plot(t_vec,rad2deg(beta_hist(3,:)))
# hold off
#
# figure
# hold on
# title('Triad Error and MEKF Error (Attitude)')
# plot(t_vec,rad2deg(angle_error),'.')
# plot(t_vec,rad2deg(ekf_angle_error))
# ylabel('Attitude Error (deg)')
# xlabel('Time (s)')
# legend('Triad','MEKF')
# hold off
# %
# % figure
# % hold on
# % plot(t_vec,ekf_angle_error)
# % hold off
#
# %% measurement
# % figure
# % hold on
# % title('Gyro Measurement')
# % plot(t_vec,meas_hist(1,:))
# % plot(t_vec,meas_hist(2,:))
# % plot(t_vec,meas_hist(3,:))
# % legend('\omega_x','\omega_y','\omega_z')
# % hold off
#
# figure
# hold on
# sgtitle('Gyro Measurement')
#
# subplot(3,1,1)
# hold on
# plot(t_vec,rad2deg(meas_hist(1,:)))
# ylabel('\omega_x (deg/s)')
# hold off
#
# subplot(3,1,2)
# hold on
# plot(t_vec,rad2deg(meas_hist(2,:)))
# ylabel('\omega_y (deg/s)')
# hold off
#
# subplot(3,1,3)
# hold on
# plot(t_vec,rad2deg(meas_hist(3,:)))
# ylabel('\omega_z (deg/s)')
# hold off
# xlabel('Time (s)')
#
#
# figure
# hold on
# title('MEKF \omega')
# plot(t_vec,mu(5,:))
# plot(t_vec,mu(6,:))
# plot(t_vec,mu(7,:))
# legend('\omega_x','\omega_y','\omega_z')
# hold off
# %
# % % get errors
# %
# % y_error = abs(state_hist(1:3,:) - meas_hist);
# ekf_error = abs(state_hist - mu);
#
# w_error = ekf_error(5:7,:);
# beta_error = ekf_error(8:10,:);


# figure
# hold on
# sgtitle('Angular Velocity Error')
#
# subplot(3,1,1)
# hold on
# plot(t_vec,rad2deg(w_error(1,:)))
# ylabel('\Delta \omega_x (deg/s)')
# hold off
#
# subplot(3,1,2)
# hold on
# plot(t_vec,rad2deg(w_error(2,:)))
# ylabel('\Delta \omega_y (deg/s)')
# hold off
#
# subplot(3,1,3)
# hold on
# plot(t_vec,rad2deg(w_error(3,:)))
# ylabel('\Delta \omega_z (deg/s)')
# hold off
# xlabel('Time (s)')
#
# figure
# hold on
# sgtitle('Gyro Bias Error')
#
# subplot(3,1,1)
# hold on
# plot(t_vec,rad2deg(beta_error(1,:)))
# ylabel('\Delta b_x (deg/s)')
# hold off
#
# subplot(3,1,2)
# hold on
# plot(t_vec,rad2deg(beta_error(2,:)))
# ylabel('\Delta b_y (deg/s)')
# hold off
#
# subplot(3,1,3)
# hold on
# plot(t_vec,rad2deg(beta_error(3,:)))
# ylabel('\Delta b_z (deg/s)')
# hold off
# xlabel('Time (s)')
# %
# % figure
# % hold on
# % title('Measuremnt error')
# % plot(t_vec,y_error(1,:))
# % plot(t_vec,y_error(2,:))
# % plot(t_vec,y_error(3,:))
# % hold off
# %
# % figure
# % hold on
# % title('EKF error')
# % plot(t_vec,ekf_error(1,:))
# % plot(t_vec,ekf_error(2,:))
# % plot(t_vec,ekf_error(3,:))
# % hold off
# %
# % figure
# % hold on
# % title('EKF bias error')
# % plot(t_vec,ekf_error(4,:))
# % plot(t_vec,ekf_error(5,:))
# % plot(t_vec,ekf_error(6,:))
# % hold off

t1 = time()
for i = 1:40000
    inv(randn(9,9))
end
@show time() - t1
