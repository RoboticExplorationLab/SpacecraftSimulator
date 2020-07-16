using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
using CSV
using DataFrames
using ProgressMeter
const SD = SatelliteDynamics
const ST = SatelliteToolbox



# load in the julia and python functions
ss_sim_path =  dirname(dirname(@__FILE__))
include(joinpath(ss_sim_path,"load_julia_functions.jl"))


function sim_driver(path_to_yaml)
    """This function runs an orbital and attitude simulation based on the
    specified configuration YAML file.

    Args:
        path_to_yaml: path from SpacecraftSimulation directory to the yaml
    """
    println("eeeee")
    # load in config from the given path
    params,initial_conditions, time_params = config(path_to_yaml)
    global params

    # timing stuff
    epc_orbital = initial_conditions.epc_orbital
    dt_orbital = time_params.dt_orbital
    dt_attitude = time_params.dt_attitude
    dt_controller = time_params.dt_controller

    # this takes care of the sample ratio
    # controller_sample_ratio = Int(dt_controller/dt_attitude)

    # time vector stuff
    tf = time_params.tf
    t_vec_orbital = time_params.t_vec_orbital
    inner_loop_t_vec = time_params.inner_loop_t_vec
    t_vec_attitude = time_params.t_vec_attitude

    # pre-allocate arrays for storage
    truth = initialize_struct(truth_state_struct,time_params,initial_conditions)

    # truth.orbital_state[1]  = initial_conditions.eci_rv_0
    # truth.attitude_state[1] = [initial_conditions.ᴺqᴮ0;
    #                              initial_conditions.ω0]

    # main loop
    t1 = time()
    @showprogress "Simulating..." for kk = 1:(length(t_vec_orbital)-1)

        # derive truth quantities
        orbital_truth_struct_update!(truth,kk,epc_orbital)

        # #----------------------ATTITUDE LOOP------------------------------
        # attitude dynamics inner loop
        for jj = 1:length(inner_loop_t_vec)-1

            # index of current step n
            index_n = (kk-1)*(length(inner_loop_t_vec)-1) + jj

            # update the derived state
            attitude_truth_struct_update!(truth,kk,index_n)

            # disturbance torques
            τ = zeros(3)

            # ------------------ control law -------------------------
            sc_mag_moment = zeros(3)

            # update the attitude with RK4
            truth.attitude_state[index_n+1] =rk4_attitude(spacecraft_eom,
            epc_orbital, truth.attitude_state[index_n], sc_mag_moment, truth.B_eci[kk], τ, dt_attitude)

            # add update for the last iter of each inner attitude loop
        end
        # --------------------------- end attitude loop -----------------------


        u_thruster = zeros(3)

        # propagate orbit one step
        truth.orbital_state[kk+1] =rk4_orbital(FODE, epc_orbital, truth.orbital_state[kk],
                                      u_thruster, dt_orbital)

        # increment the time
        epc_orbital += dt_orbital


    end

    # derive quantities for the last step
    orbital_truth_struct_update!(truth,length(t_vec_orbital),epc_orbital)
    attitude_truth_struct_update!(truth,length(t_vec_orbital),length(t_vec_attitude))

    # timing
    @show time() - t1


    return sim_output = (truth=truth, t_vec_orbital = t_vec_orbital,
                         t_vec_attitude = t_vec_attitude)
    # return sim_output = (
end

path_to_yaml = "sim/config_attitude_test.yml"
sim_output = sim_driver(path_to_yaml)


r_eci = mat_from_vec(sim_output.truth.r_eci)
B_eci = mat_from_vec(sim_output.truth.B_eci)
w = mat_from_vec(sim_output.truth.ω)
# plot(sim_output.t_vec_orbital,vec(r_eci[1,:]))

# plot orbital motion
plot(vec(r_eci[1,:])/1e3,vec(r_eci[2,:])/1e3,vec(r_eci[3,:])/1e3,title = "Orbital Motion",
label = "",xlabel = "ECI X (km)",ylabel = "ECI Y (km)", zlabel = "ECI Z (km)")

#
plot(sim_output.t_vec_orbital,vec(B_eci[1,:]))
plot!(sim_output.t_vec_orbital,vec(B_eci[2,:]))
plot!(sim_output.t_vec_orbital,vec(B_eci[3,:]))

plot(sim_output.t_vec_attitude,vec(w[1,:]))
plot!(sim_output.t_vec_attitude,vec(w[2,:]))
plot!(sim_output.t_vec_attitude,vec(w[3,:]),xlim = (0,21))
