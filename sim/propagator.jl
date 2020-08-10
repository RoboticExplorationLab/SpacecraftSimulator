using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, MATLAB

using Infiltrator
# using StaticArrays
# using CSV
# using DataFrames
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

    # load in the sim config yaml from the given path
    initial_conditions, time_params = config(path_to_yaml)

    # pre-allocate structs/arrays for storage
    truth = initialize_struct(truth_state_struct,time_params,initial_conditions)

    # main loop
    t1 = time()
    @showprogress "Simulating..." for orb_ind = 1:(length(time_params.t_vec_orbital)-1)

        # use orbital state to populate the truth struct
        orbital_truth_struct_update!(truth,orb_ind)

        # #----------------------ATTITUDE LOOP------------------------------
        # attitude dynamics inner loop
        for att_ind = 1:length(time_params.inner_loop_t_vec)-1

            # index of current step n
            index_n = (orb_ind-1)*(length(time_params.inner_loop_t_vec)-1) + att_ind

            # update the derived state
            attitude_truth_struct_update!(truth,orb_ind,index_n)

            # disturbance torques
            τ = zeros(3)

            # ------------------ control law -------------------------
            sc_mag_moment = zeros(3)

            truth.ᴺqᴮ[index_n+1], truth.ω[index_n+1] =rk4_attitude(spacecraft_eom,
            truth.epc_orbital, [truth.ᴺqᴮ[index_n];truth.ω[index_n]], sc_mag_moment, truth.B_eci[orb_ind], τ, time_params.dt_attitude)


        end
        # --------------------------- end attitude loop -----------------------

        u_thruster = zeros(3)

        # propagate orbit one step
        truth.r_eci[orb_ind+1], truth.v_eci[orb_ind+1], truth.epc_orbital = rk4_orbital(FODE, truth.epc_orbital, [truth.r_eci[orb_ind]; truth.v_eci[orb_ind]],
                                      u_thruster, time_params.dt_orbital)


    end

    # derive quantities for the last step
    orbital_truth_struct_update!(truth,length(time_params.t_vec_orbital))
    attitude_truth_struct_update!(truth,length(time_params.t_vec_orbital),length(time_params.t_vec_attitude))

    # timing
    @show time() - t1


    return sim_output = (truth=truth, t_vec_orbital = time_params.t_vec_orbital,
                         t_vec_attitude = time_params.t_vec_attitude)
    # return sim_output = (
end

path_to_yaml = "sim/config_attitude_test.yml"
sim_output = sim_driver(path_to_yaml)

ω = mat_from_vec(sim_output.truth.ω)

mat"
figure
hold on
plot($sim_output.t_vec_attitude,$ω')
"








# r_eci = mat_from_vec(sim_output.truth.r_eci)
# B_eci = mat_from_vec(sim_output.truth.B_eci)
# w = mat_from_vec(sim_output.truth.ω)
# # plot(sim_output.t_vec_orbital,vec(r_eci[1,:]))
#
# # plot orbital motion
# plot(vec(r_eci[1,:])/1e3,vec(r_eci[2,:])/1e3,vec(r_eci[3,:])/1e3,title = "Orbital Motion",
# label = "",xlabel = "ECI X (km)",ylabel = "ECI Y (km)", zlabel = "ECI Z (km)")
#
# #
# plot(sim_output.t_vec_orbital,vec(B_eci[1,:]),label = "B_x")
# plot!(sim_output.t_vec_orbital,vec(B_eci[2,:]),label = "B_y")
# plot!(sim_output.t_vec_orbital,vec(B_eci[3,:]), label = "B_z",title = "Magnetic Field in ECI (T)")
#
# plot(sim_output.t_vec_attitude,vec(w[1,:]))
# plot!(sim_output.t_vec_attitude,vec(w[2,:]))
# plot!(sim_output.t_vec_attitude,vec(w[3,:]),xlim = (0,21))
