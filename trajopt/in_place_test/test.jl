using LinearAlgebra, BenchmarkTools, MATLAB

function ode(t,x,params)
    # euler's
    return params.invJ*(-cross(x,params.J*x))
end

function ode!(t,xdot,x,params)
    # euler's inplace version
    xdot .= params.invJ*(-cross(x,params.J*x))
end


function vanila_rk4(f,t_n,x_n,dt,params)
    # standard rk4

    k1 = dt*f(t_n,x_n,params)
    k2 = dt*f(t_n+dt/2, x_n+k1/2,params)
    k3 = dt*f(t_n+dt/2, x_n+k2/2,params)
    k4 = dt*f(t_n+dt, x_n+k3,params)

    x_np1 = (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

        return x_np1
end

# function inplace_rk4(f!,t_n,xdot,x_n,dt,params)
#     # inplace version of RK4
#     x_np1 = zeros(3)
#     x_np1 .= x_n
#
#     # form k1
#     f!(t_n,xdot,x_n,params)
#     x_np1 .+= (1/6)*dt*xdot
#
#     # form k2
#     f!(t_n+dt/2,xdot, x_n+dt*xdot/2,params)
#     x_np1 .+= (1/3)*dt*xdot
#
#     # form k3
#     f!(t_n+dt/2,xdot, x_n+dt*xdot/2,params)
#     x_np1 .+= (1/3)*dt*xdot
#
#     # form k4
#     f!(t_n+dt,xdot, x_n+dt*xdot,params)
#     x_np1 .+= (1/6)*dt*xdot
#
#         return x_np1
# end

function inplace_rk4!(f!,t_n,xdot,X,i,dt,params)
    # inplace rk4

    # x_n
    X[i+1]=X[i]

    # form k1
    f!(t_n,xdot,X[i],params)
    X[i+1] += (1/6)*dt*xdot

    # form k2
    f!(t_n+dt/2,xdot, X[i]+dt*xdot/2,params)
    X[i+1] += (1/3)*dt*xdot

    # form k3
    f!(t_n+dt/2,xdot, X[i]+dt*xdot/2,params)
    X[i+1] += (1/3)*dt*xdot

    # form k4
    f!(t_n+dt,xdot, X[i]+dt*xdot,params)
    X[i+1] += (1/6)*dt*xdot

    # @infiltrate
end

function integrate()

    # inertia
    J = diagm([1;2;3])
    invJ = inv(J)
    params = (J=J,invJ = invJ)

    # initial conditions
    x0 = [.2;7;.2]

    # sample time
    dt = .01

    # final time
    tf = 10

    # time vec
    t_vec = 0:dt:tf

    # number of steps
    N = length(t_vec)

    # pre-allocate X
    X = fill(zeros(3),N)

    # set initial condition
    X[1] = x0

    # integrate
    for i = 1:N-1
        X[i+1] = vanila_rk4(ode,t_vec[i],X[i],dt,params)
    end

        return X
end

function inplace_integrate()

    # inertia
    J = diagm([1;2;3])
    invJ = inv(J)
    params = (J=J,invJ = invJ)

    # initial conditions
    x0 = [.2;7;.2]

    # sample time
    dt = .01

    # final time
    tf = 10

    # time vec
    t_vec = 0:dt:tf

    # number of steps
    N = length(t_vec)

    # pre-allocate X
    X = fill(zeros(3),N)

    # set initial condition
    X[1] .= x0

    # allocate xdot
    xdot = zeros(3)

    for i = 1:N-1
        # X[i+1] = inplace_rk4(ode!,t_vec[i],xdot,X[i],dt,params)
        inplace_rk4!(ode!,t_vec[i],xdot,X,i,dt,params)
    end

        return X
end


@btime X = integrate()

@btime X = inplace_integrate()

X = mat_from_vec(X)

mat"plot($X')"


@btime J = diagm([1;2;3])
@btime invJ = inv(J)
@btime params = (J=J,invJ = invJ)

@btime x0 = [.2;7;.2]

@btime dt = .01

@btime tf = 10

@btime t_vec = 0:dt:tf

@btime N = length(t_vec)

@btime X = fill(zeros(3),N)

@btime X[1] = x0;

for i = 1:N-1
    X[i+1] = vanila_rk4(ode,t_vec[i],X[i],dt,params)
end
