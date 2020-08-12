using LinearAlgebra, SatelliteDynamics, Test


# load in the julia and python functions
ss_sim_path =  dirname(dirname(@__FILE__))
include(joinpath(ss_sim_path,"load_julia_functions.jl"))


function SD_propagator(path_to_yaml)

    initial_conditions, time_params = config(path_to_yaml)

    # initial time
    epc0 = params.initial_conditions.epc_orbital
    # Declare initial state in terms of osculating orbital elements
    oe0 = params.initial_conditions.oe0

    # Convert osculating elements to Cartesian states
    eci0 = sOSCtoCART(oe0, use_degrees=false)

    # check everything looks good
    @assert isequal(eci0,params.initial_conditions.eci_rv_0)

    # Set the propagation end time to one orbit period after the start
    # T    = 10*orbit_period(oe0[1])
    # T = 24*3600*days
    epcf = epc0 + params.time_params.tf

    # Initialize State Vector
    orb  = SD.EarthInertialState(epc0, eci0, dt=1,
                mass=100.0, n_grav=params.grav_deg, m_grav=params.grav_order,
                drag=false, srp=false,
                moon=false, sun=false,
                relativity=false
    )

    # Propagate the orbit
    t, epc, eci = sim!(orb, epcf)

    return eci
end


path_to_yaml = "sim/config_orbital_test.yml"

eci_SD = SD_propagator(path_to_yaml)

sim_output = sim_driver(path_to_yaml)

r_eci = mat_from_vec(sim_output.truth.r_eci)
v_eci = mat_from_vec(sim_output.truth.v_eci)
eci_mine = [r_eci;v_eci]

error_mat = eci_SD - eci_mine

@testset "SD propagator reg test" begin

        @test isapprox(norm(error_mat),0.0,rtol = 1e-6)

end
