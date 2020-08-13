using LinearAlgebra, JLD2

@load "B.jld2" B_mat B_save t_vec

# first let's do interpolation
function interp1(t,B_save,input_t)
    # sample rate of B_save vec
    Δt = t[2]-t[1]

    # get the lower and upper bounds
    lower_idx = Int(floor(input_t/Δt))
    # upper_idx = Int(ceil(input_t+1e-10/Δt))

    # how far are we between the bounds ∈[0,1]
    δt = ((input_t) - Δt*lower_idx)/Δt

    # lower and upper function values
    f_lower = B_save[lower_idx+1]
    f_upper = B_save[lower_idx+2]

        return (f_lower + δt*(f_upper - f_lower))
end


J = [1.959e-4 2016.333e-9 269.176e-9;
           2016.333e-9 1.999e-4 2318.659e-9;
           269.176e-9 2318.659e-9 1.064e-4]

invJ = inv(J)

alpha0 = 1.0

max_moments = [8.8e-3;1.373e-2;8.2e-3]
dJ_tol = 1e-2
params = (alpha0 = alpha0, t_vec = t_vec, B_save = B_save, J = J, invJ = invJ,
max_moments = max_moments,dJ_tol = dJ_tol)

# x0 = [deg2rad(120)*normalize(randn(3));zeros(3)]
x0 = vec([-1.4505896726335037, 0.5147044464043511, -1.420337910297834, 0.0, 0.0, 0.0])

xg = zeros(6)

N = 120

# N = 60;
u0 = 1*zeros(3,N-1);

Q = 1*eye(6)
Q[4:6,4:6] = .01*Q[4:6,4:6]
R = 1*eye(3)
Qf = 100*eye(6)
Qf[4:6,4:6] = 10*Qf[4:6,4:6]
dt = 2.5

xhist, uhist, K = iLQRsimple_B(x0, xg, u0, Q, R, Qf, dt,params)


mat"
figure
hold on
plot($xhist(1:3,:)')
"

uplot = diagm(params.max_moments)*tanh.(uhist)
mat"
figure
hold on
plot($uplot')
"
