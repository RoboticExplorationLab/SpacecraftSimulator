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
