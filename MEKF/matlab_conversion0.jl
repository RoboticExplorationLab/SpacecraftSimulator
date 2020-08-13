using LinearAlgebra

mutable struct EKF_struct
    mu::Array{Array{Float64,1},1}
    Sigma::Array{Array{Float64,2},1}
end

function dynODE(t_n,x_n,tau,params)

    # J = params.J
    # invJ = params.invJ

    # x_n = x_n[:];
    q = x_n[1:4]
    w = x_n[5:7]
    # beta = x_n(4:6);

    q_dot = .5*qdot(q,[w;0])

    wdot = params.invJ*(tau - cross(w,params.J*w))

    beta_dot = zeros(3)

    Xdot = [q_dot;
            wdot;
            beta_dot]

    return Xdot

end

function g(x,r1n, r2n)

    # x = x[:];
    q = x[1:4]
    w = x[5:7]

    beta = x[8:10]

    gyro = w+beta

    n_Q_b = dcm_from_q(q)
    # b_Q_n = n_Q_b'

    r1b = transpose(n_Q_b)*r1n
    r2b = transpose(n_Q_b)*r2n

    y = [gyro;r1b;r2b]

    return y

end

function add_noise_to_y(y,sqrtR)
# add measurement noise t
# y = y[:];
gyro = y[1:3]
r1b = y[4:6]
r2b = y[7:9]


noisy_gyro = params.sensors.offsets.gyroscope_offset*gyro + fast_mvnrnd(sqrtR[1:3,1:3])

# sun sensor
noisy_r1b = params.sensors.offsets.sun_sensor_offset*S03_noise(r1b,sqrtR[4:6,4:6])

# magnetometer
noisy_r2b = params.sensors.offsets.magnetometer_offset*S03_noise(r2b,sqrtR[7:9,7:9])

noisy_y = [noisy_gyro;
           noisy_r1b;
           noisy_r2b]

return noisy_y

end


function S03_noise(vector, sqrt_covariance)
    # here we add noise S03 style

    # noise axis angle vector
    phi_noise = fast_mvnrnd(sqrt_covariance)

    # apply the noise rotation to the vector
    noisy_vector = skew_expm(hat(phi_noise))*vector

    return noisy_vector
end

function slim_ekf_predict(mu,Sigma,tau,w_gyro,params,dt,Q)

    # unpack state
    # mu = mu[:];
    q = mu[1:4]
    beta = mu[5:7]

    # predict
    q_next = q ⊙ q_from_phi((w_gyro - beta)*dt)
    beta_next = beta

    mu_predict = [q_next;
                  beta_next]

    # jacobians
    dphi_phi = I- hat((w_gyro - beta)*dt)
    dphi_dbeta = -.5*dt*eye(3)
    dbeta_beta = eye(3)

    # dynamics jacobian (d_f/d_state)
    A = [dphi_phi      dphi_dbeta;
         zeros(3,3)    dbeta_beta]

    # covariance prediction
    # @infiltrate
    # error()
    Sigma_predict = A*Sigma*A' + Q

    return mu_predict, Sigma_predict
end

function slim_ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R)

    """yt is just r1b and r2n"""

    # predicted quaternion
    q_predict = mu_predict[1:4]

    # dcm
    n_Q_bt = dcm_from_q(q_predict)

    # use the predicted attitude to generate predicted measurements
    r1b_t = transpose(n_Q_bt)*r1n
    r2b_t = transpose(n_Q_bt)*r2n


    # predicted measurment
    yhat = [r1b_t;r2b_t]

    # measurement difference
    z = yt - yhat

    # measurement Jacobian
    C = [2*hat(r1b_t) zeros(3,3) ;
         2*hat(r2b_t) zeros(3,3) ]

    # covariance of the innovation
    S = C*Sigma_predict*transpose(C) + R

    # kalman gain
    # L = Sigma_predict*C'*inv(S)
    L = Sigma_predict*transpose(C)/S

    # innovation
    delta_mu = L*z

    # update quaternion multiplicatively
    q_update = q_predict ⊙ q_from_g(delta_mu[1:3])

    # update bias additively
    beta_update = mu_predict[5:7] + delta_mu[4:6]

    # updated mu
    mu_update = [ q_update ;
                  beta_update]

    # Sigma update
    Sigma_update = (eye(6) - L*C)*Sigma_predict

    return mu_update, Sigma_update
end


function ekf_predict(mu,Sigma,tau,params,dt,Q)

    # tn = 0;

    # unpack state
    # mu = mu[:];
    q = mu[1:4]
    w = mu[5:7]
    beta = mu[8:10]

    # predict
    q_next = q ⊙ q_from_phi(w*dt)
    wdot = params.sensors.Jinv_sample*(tau-cross(w,params.sensors.J_sample*w))
    w_next = w + wdot*dt
    beta_next = beta

    mu_predict = [q_next;
                  w_next;
                  beta_next]

    # jacobians
    dphi_phi = I- hat(w*dt)
    dphi_dw = .5*dt*eye(3)
    dw_w = (I + dt*params.sensors.Jinv_sample*(hat(params.sensors.J_sample*w) - hat(w)*params.sensors.J_sample))
    dbeta_beta = eye(3)

    # dynamics jacobian (d_f/d_state)
    A = [dphi_phi        dphi_dw          zeros(3,3);
         zeros(3,3)      dw_w                zeros(3,3);
         zeros(3,3)      zeros(3,3)          dbeta_beta]

    # covariance prediction
    Sigma_predict = A*Sigma*A' + Q

    return mu_predict, Sigma_predict
end


function ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R)

    # predicted quaternion
    q_predict = mu_predict[1:4]

    # dcm
    n_Q_bt = dcm_from_q(q_predict)

    # use the predicted attitude to generate predicted measurements
    r1b_t = transpose(n_Q_bt)*r1n
    r2b_t = transpose(n_Q_bt)*r2n


    # predicted measurment
    yhat = g(mu_predict,r1n, r2n)

    # measurement difference
    z = yt - yhat

    # measurement Jacobian
    C = [zeros(3,3)     I          I  ;
         2*hat(r1b_t) zeros(3,3) zeros(3,3);
         2*hat(r2b_t) zeros(3,3) zeros(3,3)]

    # covariance of the innovation
    S = C*Sigma_predict*transpose(C) + R

    # kalman gain
    # L = Sigma_predict*C'*inv(S)
    L = Sigma_predict*transpose(C)/S

    # innovation
    delta_mu = L*z

    # update quaternion multiplicatively
    q_update = q_predict ⊙ q_from_g(delta_mu[1:3])

    # update angular velocity and bias additively
    w_update = mu_predict[5:7] + delta_mu[4:6]
    beta_update = mu_predict[8:10] + delta_mu[7:9]

    # updated mu
    mu_update = [ q_update ;
                  w_update;
                  beta_update]

    # Sigma update
    Sigma_update = (eye(9) - L*C)*Sigma_predict

    return mu_update, Sigma_update
end

function MEKF(mu,Sigma,tau,yt,Q,R,dt,r1n,r2n,params)

# predict
mu_predict,Sigma_predict = ekf_predict(mu,Sigma,tau,params,dt,Q)

# innovate/update
mu_update,Sigma_update = ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R)

    return mu_update, Sigma_update
end

function slim_MEKF(mu,Sigma,tau,w_gyro,yt,Q,R,dt,r1n,r2n,params)

# predict
mu_predict,Sigma_predict = slim_ekf_predict(mu,Sigma,tau,w_gyro,params,dt,Q)

# innovate/update
mu_update,Sigma_update = slim_ekf_innovate(mu_predict,Sigma_predict,yt,r1n,r2n,R)

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

    k1 = h*ODE(tn,xn,u,sc)
    k2 = h*ODE(tn + h/2,xn + k1/2,u,sc)
    k3 = h*ODE(tn + h/2,xn + k2/2,u,sc)
    k4 = h*ODE(tn + h,xn + k3,u,sc)

    y_np1 = xn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    return y_np1


end

function sample_inertia(J,deg_std,scale_std)

    # take eigen decomposition
    eigen_decomp = eigen(J)
    D = diagm(eigen_decomp.values)
    S = eigen_decomp.vectors

    # scale the moments
    D_sample = D*diagm(ones(3)+ scale_std*randn(3))

    # rotate them
    R_sample = exp(hat(deg2rad(deg_std)*randn(3)))
    J_sample = R_sample*S*D_sample*S'*R_sample'

    return J_sample
end


# %% Sim
function run_MEKF()
# % inertia
J = [   1.959e-4    2016.333e-9     269.176e-9;
     2016.333e-9       1.999e-4    2318.659e-9;
      269.176e-9    2318.659e-9       1.064e-4]

J_sample = sample_inertia(J,10.0,.1)

sample_rate = 10

sun_sensor = (offset_std = deg2rad(3),noise_std = deg2rad(10))
magnetometer = (offset_std = deg2rad(3),noise_std = deg2rad(10))
gyroscope = (offset_std = deg2rad(3),bias_walk_std = 1e-8,noise_std = deg2rad(sqrt(.07^2*sample_rate)))

# generate offsets
sun_sensor_offset = skew_expm(hat(sun_sensor.offset_std*randn(3)))
magnetometer_offset = skew_expm(hat(magnetometer.offset_std*randn(3)))
gyroscope_offset = skew_expm(hat(gyroscope.offset_std*randn(3)))

# store offsets
offsets = (sun_sensor_offset = sun_sensor_offset,
           magnetometer_offset = magnetometer_offset,
           gyroscope_offset = gyroscope_offset)

# store all sensor information here
sensors = (sun_sensor = sun_sensor,
           magnetometer = magnetometer,
           gyroscope = gyroscope,
           offsets = offsets,
           J_sample = J_sample,
           Jinv_sample = inv(J_sample))

# sc = (sensors = sensors, J_sampled = )

# params.J = J; params.invJ = invJ;
global params = (J=J,invJ = inv(J), sensors = sensors )

# % initial conditions
q0 = randq()
w0 = deg2rad(3)*normalize(randn(3))
beta0 = deg2rad(3)*normalize(randn(3))

# % two known vectors expressed in ECI
r1n = normalize(randn(3));
r2n = normalize(randn(3));

# time stuff
dt = .1
tf = 400.0
t_vec = 0:dt:tf

# pre-allocate
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
slim_ekf_angle_error = zeros(length(t_vec));

# % Process noise covariance
Q = zeros(10,10);
Q[8:10,8:10] = 1e-9*eye(3);
Q_ekf = .000001*eye(9);
Q_ekf[7:9,7:9] = Q[8:10,8:10];

sqrtQ = sqrt(Q)
# % Sensor noise covariance
# R_gyro = (deg2rad(1))^2*eye(3);
# R_sun_sensor = (deg2rad(5))^2*eye(3);
# R_magnetometer = (deg2rad(5))^2*eye(3);
R_gyro = params.sensors.gyroscope.noise_std^2*eye(3)
R_sun_sensor = params.sensors.sun_sensor.noise_std^2*eye(3)
R_magnetometer = params.sensors.magnetometer.noise_std^2*eye(3)
R = zeros(9,9)
R[1:3,1:3] = R_gyro
R[4:6,4:6] = R_sun_sensor
R[7:9,7:9] = R_magnetometer
R_ekf = 4*R
sqrtR = sqrt(R)

# % external torque on spacecraft
tau = zeros(3)

mekf_time = 0.0

# % MEKF initialize
EKF = EKF_struct(copy(fill(zeros(10),length(t_vec))),copy(fill(zeros(9,9),length(t_vec))))
slim_EKF = EKF_struct(copy(fill(zeros(7),length(t_vec))),copy(fill(zeros(6,6),length(t_vec))))
# mu = state_hist;
# Sigma = cell(length(t_vec),1);
EKF.mu[1] = [state_hist[1:7,1];
             zeros(3)]
EKF.Sigma[1] = .1*eye(9);

slim_EKF.mu[1] = [state_hist[1:4,1];
                  zeros(3)]
slim_EKF.Sigma[1] = .1*eye(6)

slim = (Q = diagm([1e-6*ones(3);1e-9*ones(3)]),R = 1.5*diagm([diag(R_sun_sensor);diag(R_magnetometer)]))

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

    # triad angle error
    angle_error[kk+1] = q_angle_error(state_hist[1:4,kk+1],n_q_b_est);

    # MEKF
    EKF.mu[kk+1], EKF.Sigma[kk+1] = MEKF(EKF.mu[kk],EKF.Sigma[kk],tau,meas_hist[:,kk+1],Q_ekf,R_ekf,dt,r1n,r2n,params)
    slim_EKF.mu[kk+1], slim_EKF.Sigma[kk+1] = slim_MEKF(slim_EKF.mu[kk],slim_EKF.Sigma[kk],tau,meas_hist[1:3,kk],meas_hist[4:9,kk+1],slim.Q,slim.R,dt,r1n,r2n,params)

    # % MEKF angle error
    ekf_angle_error[kk+1] = q_angle_error(state_hist[1:4,kk+1],EKF.mu[kk+1][1:4]);
    slim_ekf_angle_error[kk+1] = q_angle_error(state_hist[1:4,kk+1],slim_EKF.mu[kk+1][1:4]);
end

mat"
figure
hold on
plot($t_vec,rad2deg($ekf_angle_error))
plot($t_vec,rad2deg($slim_ekf_angle_error))
plot($t_vec,rad2deg($angle_error),'.')
legend('MEKF','Slim MEKF','Triad')
xlabel('Time (s)')
ylabel('Pointing Error (deg)')
title('Triad and MEKF Errors')
"

# @show mekf_time

return t_vec, meas_hist, state_hist, EKF, slim_EKF
end


t_vec, meas_hist, state_hist, EKF, slim_EKF = run_MEKF()


ekf_mu= mat_from_vec(EKF.mu)
slim_ekf_mu= mat_from_vec(slim_EKF.mu)
# meas_hist = mat_from_vec(meas_hist)

w_hist = state_hist[5:7,:]
b_hist = state_hist[8:10,:]
ekf_w = ekf_mu[5:7,:]
ekf_b = ekf_mu[8:10,:]
gyro_w = meas_hist[1:3,:]
slim_ekf_w = gyro_w - slim_ekf_mu[5:7,:]
slim_ekf_b = slim_ekf_mu[5:7,:]

N = size(w_hist,2)
ekf_w_error = zeros(N)
ekf_b_error = zeros(N)
slim_ekf_w_error = zeros(N)
slim_ekf_b_error = zeros(N)
for i = 1:N

    ekf_w_error[i] = norm(ekf_w[:,i] - w_hist[:,i])
    slim_ekf_w_error[i] = norm(slim_ekf_w[:,i] - w_hist[:,i])
    ekf_b_error[i] = norm(ekf_b[:,i] - b_hist[:,i])
    slim_ekf_b_error[i] = norm(slim_ekf_b[:,i] - b_hist[:,i])

end

mat"
figure
hold on
title('Angular Velocity Error')
plot($t_vec,rad2deg($ekf_w_error),'.')
plot($t_vec,rad2deg($slim_ekf_w_error),'.')
legend('MEKF','slim MEKF')
ylabel('Angular Velocity Error (deg/s)')
xlabel('Time (s)')
xlim([10,400])
hold off
"

mat"
figure
hold on
title('Gyro Bias Error')
plot($t_vec,rad2deg($ekf_b_error))
plot($t_vec,rad2deg($slim_ekf_b_error))
legend('MEKF','slim MEKF')
xlim([10,400])
hold off
"

# mat"figure
# hold on
# title('Angular Velocity')
# % plot($t_vec,rad2deg($slim_ekf_w'))
# plot($t_vec,rad2deg($ekf_w'))
# plot($t_vec,rad2deg($w_hist'))
# hold off
# "
# mat"figure
# hold on
# title('Gyro Bias')
# % plot($t_vec,rad2deg($slim_ekf_w'))
# plot($t_vec,rad2deg($ekf_b'))
# plot($t_vec,rad2deg($b_hist'))
# hold off
# "
