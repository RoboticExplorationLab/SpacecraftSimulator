
using LinearAlgebra, BenchmarkTools

f_x(x) = f(t,x,u,params)

function f(t,x,u,params)

    # x is dim2
    # u is dim2

    J = params.J


    A = u*[1 3; 4.0 5]

    return A*x
end

# function A_fun(t,x,u,params)
#     return A_fx(x)
# end


function test_me()



    u = 1.0
    params = (J = diagm([1;2;3]), aa = 4.3)

    t = 32.0

    A_fx = x-> ForwardDiff.jacobian(x -> f(t,x,u,params),x)


    x = randn(2)

    @show ForwardDiff.jacobian(x -> f(t,x,u,params),x)

    @show A_fx(x)

    u = 2.0

    @time A_fx(x)

    # @show A_fun(t,x,3.0,params)
end


test_me()
