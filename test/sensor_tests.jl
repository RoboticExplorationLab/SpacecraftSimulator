using Test

# load in the julia and python functions
ss_sim_path =  dirname(dirname(@__FILE__))
include(joinpath(ss_sim_path,"load_julia_functions.jl"))

path_to_yaml =  "sim/config_attitude_test.yml"
params,initial_conditions, time_params = config(path_to_yaml)
global params

@testset "sun_body" begin
    let
        for i = 1:1000

            # generate time
            epc = params.initial_conditions.epc_orbital

            # generate attitude
            q = randq()
            ᴺQᴮ = dcm_from_q(q)

            # spacecraft position vector
            r_eci = (R_EARTH+550000)*normalize(randn(3))

            # sun position
            r_sun_eci = sun_position(epc + 100000*randn())

            # check eclipse
            eclipse = eclipse_check(r_eci, r_sun_eci)

            # sun vector in the body frame
            sun_body = transpose(ᴺQᴮ)*normalize(r_sun_eci - r_eci)

            sun_body2 = sun_body_normalized(r_sun_eci,r_eci,ᴺQᴮ)

            I_flux = sun_flux(r_sun_eci,r_eci,ᴺQᴮ,eclipse)

            # @show I_flux

            sun_body3 = s_body_from_I(I_flux)

            if !eclipse
                @test isapprox(sun_body,sun_body2,rtol = 1e-6)
                @test isapprox(sun_body,sun_body3,rtol = 1e-6)
            else
                @test isapprox(I_flux,zeros(6),rtol = 1e-6)
            end
            

        end
    end
end
