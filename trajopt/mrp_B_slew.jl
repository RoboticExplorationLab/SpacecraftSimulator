using LinearAlgebra, JLD2, MATLAB


@load "B.jld2" B_mat B_save t_vec

include(joinpath(dirname(@__FILE__),"sc_B_dynamics.jl"))
include(joinpath(dirname(@__FILE__),"ilqrsolver.jl"))
# include(joinpath(@__FILE_))



J = convert(Array{Float32,2},[1.959e-4 2016.333e-9 269.176e-9;
           2016.333e-9 1.999e-4 2318.659e-9;
           269.176e-9 2318.659e-9 1.064e-4])

invJ = inv(J)

alpha0 = 1

max_moments = convert(Array{Float32,1},[8.8e-3;1.373e-2;8.2e-3])
dJ_tol = 1e-3
params = (alpha0 = alpha0, t_vec = t_vec, B_save = convert(Array{Array{Float32,1},1},B_save), J = J, invJ = invJ,
max_moments = max_moments,dJ_tol = dJ_tol)

function runit()
x0 = [p_from_phi(deg2rad(160)*normalize(randn(3)));zeros(3)]
# x0 = vec([-0.569111, 0.08157139, -0.05284032, 0.0, 0.0, 0.0])
# x0 = vec([-1.4505896726335037, 0.5147044464043511, -1.420337910297834, 0.0, 0.0, 0.0])
# x0 = vec([-0.17620427335985153, -0.2757563693297543, -0.4756509352005231, 0.0, 0.0, 0.0])
# x0 = vec([0.2959220196740442, -0.06054147908212846, 0.4920347761245298, 0.0, 0.0, 0.0])
# x0 = vec([0.53753614, -0.1159876, 0.17588383, 0.0, 0.0, 0.0])
xg = zeros(6)

N = 25

# N = 60;
u0 = 1*zeros(3,N-1);

Q = 1*eye(6)
Q[4:6,4:6] = .01*Q[4:6,4:6]
R = .1*eye(3)
Qf = 100*eye(6)
Qf[4:6,4:6] = 10*Qf[4:6,4:6]
dt = 12
tf = N*dt


# convert everything to Float32
x0 = convert(Array{Float32,1},x0)
xg = convert(Array{Float32,1},xg)
u0 = convert(Array{Float32,2},u0)
Q = convert(Array{Float32,2},Q)
R = convert(Array{Float32,2},R)
Qf = convert(Array{Float32,2},Qf)


dt = convert(Float32,dt)

# B_eci = interp1(params.t_vec,params.B_save,11.0)

# xhist, uhist, Khist = iLQRsimple_B(x0, xg, u0, Q, R, Qf, dt,params)
xhist, uhist, Khist = iLQRsimple_B3(x0, xg, Q, R, Qf, N,dt,params)
println("solved")

# Khist = mat_from_vec(Khist)
uhist = mat_from_vec(uhist)
xhist = mat_from_vec(xhist)
angle_error = zeros(size(xhist,2))
ω_norm = zeros(size(xhist,2))
for i = 1:size(xhist,2)
    angle_error[i] = rad2deg(norm(phi_from_p(xhist[1:3,i])))
    ω_norm[i] = rad2deg(norm(xhist[4:6,i]))
end

t_plot = 0:dt:(N-1)*dt


# mat"
# figure
# hold on
# title('Pointing Error')
# ylabel('Degrees')
# xlabel('Time (s)')
# plot($t_plot,$angle_error)
# hold off
# "
#
# mat"
# figure
# hold on
# title('Angular Velocity')
# ylabel('Degrees/second')
# xlabel('Time (s)')
# plot($t_plot,$ω_norm)
# hold off
# "
#
#
#
# mat"
# figure
# hold on
# title('MRP')
# plot($t_plot,$xhist(1:3,:)')
# "
#
# uplot = diagm(params.max_moments)*tanh.(uhist)
# mat"
# figure
# hold on
# title('Magnetic Moment')
# stairs($t_plot(1:end-1),$uplot')
# "

# Khist = vec_from_mat(Khist)
uhist = vec_from_mat(uhist)
xhist = vec_from_mat(xhist)

x_t = t_plot
u_t = t_plot[1:end-1]


# now lets try a sim
function controller(u_t,uhist,x_t,xhist,Khist,t,x)

    if t>u_t[end]
        uplan = zeros(3)
    else
        uplan = interp1(u_t,uhist,t)
    end

    if t>x_t[end]
        xplan = zeros(6)
        K     = Khist[end]
    else
        xplan = interp1(x_t,xhist,t)
        K     = interp1(u_t,Khist,t)
    end

    # if
    # K     = interp1(u_t,Khist,t)

    # @infiltrate
    # error()
    dx = (x - xplan)

    u = uplan - K*dx

    return u
end


dt = 1.0
N = 300
t_vec = 0:dt:(N)*dt

X = fill(zeros(6),N)

X[1] = x0

process_noise = sqrt(diagm([.000001*ones(3);.0000000001*ones(3)]))

point_error = zeros(N-1)

for i = 1:N-1
    t = t_vec[i]

    u = controller(u_t,uhist,x_t,xhist,Khist,t,X[i]+10*process_noise*randn(6))

    # add noise on control input
    u = skew_expm(hat(deg2rad(20)*normalize(randn(3))))*u

    # add noise on state
    X[i]+= process_noise*randn(6)

    X[i+1] = rk4(sc_b_dynamics_xdot,t,X[i],u,params,dt)

    point_error[i] = rad2deg(norm(phi_from_p(X[i][1:3])))

end

X = mat_from_vec(X)

mat"
figure
hold on
plot($point_error')
hold off
"
end

runit()
