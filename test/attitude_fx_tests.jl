using LinearAlgebra, Test

@testset "randq" begin
    let
        for i = 1:100
            q = randq()

            @test isequal(ndims(q),1)
            @test isequal(size(q,1),4)
            @test isequal(length(q),4)

            @test isapprox(norm(q),1.0,rtol = 1e-6)
        end
    end
end

@testset "qshorter" begin
    let
        for i = 1:100
            q = randq()
            qs = q_shorter(q)

            # rotation between the two
            b1_q_b2 = qconj(q) ⊙ qs

            # axis angle between the two
            ϕ = phi_from_q(b1_q_b2)

            @test isapprox(norm(ϕ),0.0,rtol = 1e-6)


        end
    end
end

@testset "q <-> DCM" begin
    let
        for i = 1:1000

            q = randq()
            Q = dcm_from_q(q)
            q2 = q_from_dcm(Q)

            q = q_shorter(q)
            q2 = q_shorter(q2)

            @test isapprox(q,q2,rtol = 1e-6)


        end
    end
end

@testset "pure quaternion conversions" begin
    let
        for i = 1:1000

            x = randn(3)

            x_quaternion = H()*x

            @test isapprox(x,x_quaternion[1:3],rtol = 1e-6)

            @test isapprox(x_quaternion[1:3],H()'*x_quaternion,rtol = 1e-6)

        end
    end
end


@testset "axis angle rotations" begin
    let
        for i = 1:1000

            ϕ = randn(3)

            x_old = randn(3)

            q = q_from_phi(ϕ)
            Q = dcm_from_phi(ϕ)

            x_new_expm = exp(skew_from_vec(ϕ))*x_old

            x_new_dcm = Q*x_old

            x_new_quaternion = H()'*(q ⊙ (H()*x_old) ⊙ qconj(q))

            @test isapprox(x_new_expm,x_new_dcm,rtol = 1e-6)

            @test isapprox(x_new_expm,x_new_quaternion,rtol = 1e-6)

        end
    end
end

@testset "shorter quaternion" begin
    let
        for i = 1:1000
            q = randq()
            q = q_shorter(q)

            phi_shorter = phi_from_q(q)
            phi_longer = phi_from_q(-q)

            @test norm(phi_shorter) < norm(phi_longer)
        end
    end
end

@testset "q - p - g" begin
    let
        for i = 1:1000

            q = randq()
            g = g_from_q(q)
            p = p_from_q(q)

            q_g = q_from_g(g)
            q_p = q_from_p(p)

            @test isapprox(q_shorter(q),q_shorter(q_g),rtol = 1e-6)
            @test isapprox(q_shorter(q),q_shorter(q_p),rtol = 1e-6)

        end
    end
end


@testset "DCM conversions" begin
    let
        for i = 1:1000

            x = randn(3)

            q = randq()
            g = g_from_q(q)
            p = p_from_q(q)

            dcm_q = dcm_from_q(q)
            dcm_g = dcm_from_g(g)
            dcm_p = dcm_from_p(p)

            @test matrix_isapprox(dcm_q,dcm_g,1e-7)
            @test matrix_isapprox(dcm_q,dcm_p,1e-7)

            @test isapprox(dcm_q*x,dcm_p*x,rtol = 1e-6)
            @test isapprox(dcm_q*x,dcm_g*x,rtol = 1e-6)

        end
    end
end
