using Test


J = [1.959e-4 2016.333e-9 269.176e-9;
           2016.333e-9 1.999e-4 2318.659e-9;
           269.176e-9 2318.659e-9 1.064e-4]

invJ = inv(J)

alpha0 = .25

max_moments = [8.8e-3;1.373e-2;8.2e-3]
dJ_tol = 1e-2
global params = (alpha0 = alpha0, t_vec = t_vec, B_save = B_save, J = J, invJ = invJ,
max_moments = max_moments,dJ_tol = dJ_tol)




# here we forward diff the dynamics
@testset "dynamics jacobians" begin
    let
        for i = 1:100
            t = 403.4
            u = randn(3)

            x_test = randn(6)


            function fwd_f(x)
                xdot,A,B = sc_b_dynamics(t,x,u,params)
                return xdot
            end

            J = x -> ForwardDiff.jacobian(fwd_f,x)


            xdot,A,B = sc_b_dynamics(t,x_test,u,params)

            Afd = J(x_test)


            function fwd_u(u)
                xdot,A,B = sc_b_dynamics(t,x_test,u,params)
                return xdot
            end


            J_2 = u -> ForwardDiff.jacobian(fwd_u,u)

            Bfd = J_2(u)

            @test norm(A-Afd)<1e-13
            @test norm(B-Bfd)<1e-13
        end
    end
end
