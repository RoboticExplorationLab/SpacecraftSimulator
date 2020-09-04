ss_sim_path =  dirname(dirname(@__FILE__))

cd(joinpath(ss_sim_path,"virtual_env"))
Pkg.activate(".")


using LinearAlgebra, SatelliteDynamics, MATLAB, ProgressMeter
using Infiltrator

# keep namespaces consistent
const SD = SatelliteDynamics



# load in the julia and python functions
ss_sim_path =  dirname(dirname(@__FILE__))
include(joinpath(ss_sim_path,"load_julia_functions.jl"))


function sim_driver(path_to_yaml)
    """This function runs an orbital and attitude simulation based on the
    specified configuration YAML file.

    Args:
        path_to_yaml: path from SpacecraftSimulation directory to the yaml

    Returns:
        sim_output: named tuple
    """

    # load in the sim config yaml from the given path
    initial_conditions, time_params = config(path_to_yaml)

    # pre-allocate structs/arrays for storage
    truth = initialize_struct(truth_state_struct,time_params,initial_conditions)
    sensors = initialize_sensors_struct(time_params)
    MEKF = initialize_mekf_struct(time_params)

    # MEKF
    initialize_mekf!(MEKF,truth)

    # triad
    q_triad = fill(zeros(4),length(time_params.t_vec_attitude))

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

            # generate measurements
            sensors_update!(sensors, truth,orb_ind,index_n)

            # disturbance torques
            τ = zeros(3)

            # --------------------MEKF--------------------------------


            if index_n>1
                # mekf!(MEKF,τ,sensors,orb_ind,index_n)
                yt = [sensors.ω[index_n];
                      sensors.sun_body[index_n];
                      sensors.B_body[index_n]]

                MEKF.mu[index_n], MEKF.Sigma[index_n] = mekf(MEKF.mu[index_n-1],
                    MEKF.Sigma[index_n-1],τ,normalize(sensors.sun_eci[orb_ind]),
                    normalize(sensors.B_eci[orb_ind]),yt,MEKF.Q,MEKF.R,
                    params.time_params.dt_attitude,sensors.J,sensors.invJ)
            end
            update_mekf_struct!(MEKF,index_n)

            # -------------------TRIAD--------------------------------
            q_triad[index_n] = triad_equal_error(sensors.sun_eci[orb_ind],sensors.B_eci[orb_ind],
                                                 sensors.sun_body[index_n],sensors.B_body[index_n])
            # ------------------ control law -------------------------
            sc_mag_moment = zeros(3)

            # -------------------integration--------------------------
            truth.ᴺqᴮ[index_n+1], truth.ω[index_n+1], truth.β[index_n+1] =rk4_attitude(spacecraft_eom,
            truth.epc_orbital, [truth.ᴺqᴮ[index_n];truth.ω[index_n]],
            truth.β[index_n],sc_mag_moment, truth.B_eci[orb_ind], τ,
            time_params.dt_attitude)


        end
        # --------------------------- end attitude loop -----------------------

        u_thruster = zeros(3)

        # propagate orbit one step
        truth.r_eci[orb_ind+1], truth.v_eci[orb_ind+1], truth.epc_orbital = rk4_orbital(FODE, truth.epc_orbital, [truth.r_eci[orb_ind]; truth.v_eci[orb_ind]],
                                      u_thruster, time_params.dt_orbital)


    end

    # derive quantities for the last step
    index_n = length(time_params.t_vec_attitude)
    orb_ind = length(time_params.t_vec_orbital)
    orbital_truth_struct_update!(truth,length(time_params.t_vec_orbital))
    attitude_truth_struct_update!(truth,length(time_params.t_vec_orbital),
                                  length(time_params.t_vec_attitude))
    sensors_update!(sensors, truth,length(time_params.t_vec_orbital),
                                  length(time_params.t_vec_attitude))
    # yt = [sensors.ω[index_n];sensors.sun_body[index_n];sensors.B_body[index_n]]
    # MEKF.mu[index_n], MEKF.Sigma[index_n] = mekf(MEKF.mu[index_n-1],
    #         MEKF.Sigma[index_n-1],τ,normalize(sensors.sun_eci[orb_ind]),
    #         normalize(sensors.B_eci[orb_ind]),yt,MEKF.Q,MEKF.R,
    #         params.time_params.dt_attitude,sensors.J,sensors.invJ)

    return sim_output = (truth=truth, t_vec_orbital = time_params.t_vec_orbital,
                         t_vec_attitude = time_params.t_vec_attitude, sensors,
                         MEKF = MEKF, q_triad)

end

path_to_yaml = "sim/config_attitude_test0.yml"
sim_output = sim_driver(path_to_yaml)

q = mat_from_vec(sim_output.truth.ᴺqᴮ)
ω = mat_from_vec(sim_output.truth.ω)

q_triad = mat_from_vec(sim_output.q_triad)

β = mat_from_vec(sim_output.truth.β)
β_ekf = mat_from_vec(sim_output.MEKF.β)
mu = mat_from_vec(sim_output.MEKF.mu)
mu_w = mu[5:7,:]
gyro = mat_from_vec(sim_output.sensors.ω)


N = size(mu_w,2)
mekf_point_error = zeros(N)
triad_point_error = zeros(N)
for i = 1:N
    mekf_point_error[i] = q_angle_error(sim_output.truth.ᴺqᴮ[i],mu[1:4,i])
    triad_point_error[i] = q_angle_error(sim_output.truth.ᴺqᴮ[i],q_triad[:,i])
end

mat"
figure
hold on
plot($sim_output.t_vec_attitude,$β_ekf')
plot($sim_output.t_vec_attitude,$β')
"

mat"
figure
hold on
plot($sim_output.t_vec_attitude, rad2deg($triad_point_error),'.')
plot($sim_output.t_vec_attitude, rad2deg($mekf_point_error))
hold off
"

# mat"
# figure
# hold on
# plot($sim_output.t_vec_attitude,$q')
# plot($sim_output.t_vec_attitude,$q_triad')
# %plot($sim_output.t_vec_attitude,$mu(1:4,:)')
# "
#
#
# mat"
# figure
# hold on
# plot($sim_output.t_vec_attitude, $mekf_point_error)
# hold off
# "
mat"
figure
hold on
plot($sim_output.t_vec_attitude,$ω')
plot($sim_output.t_vec_attitude,$gyro')
"
mat"
figure
hold on
plot($sim_output.t_vec_attitude,$ω')
plot($sim_output.t_vec_attitude,$mu_w')
"
