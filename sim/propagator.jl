using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
using CSV
using DataFrames
using ProgressMeter
const SD = SatelliteDynamics
const ST = SatelliteToolbox


mutable struct truth_struct
    orbital_state   :: Array{Float64,2}
    attitude_state  :: Array{Float64,2}
    B_eci           :: Array{Float64,2}
    eclipse_hist    :: Array{Bool,1}
    r_sun_eci       :: Array{Float64,2}
    I_sun_flux      :: Array{Float64,2}
end

mutable struct sense_struct
    orbital_state   :: Array{Float64,2}
    attitude_state  :: Array{Float64,2}
    B_eci           :: Array{Float64,2}
    B_body          :: Array{Float64,2}
    eclipse_hist    :: Array{Bool,1}
end


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



    # @infiltrate
    # error()
    truth = truth_struct(zeros(6,length(t_vec_orbital)),
                         zeros(7,length(t_vec_attitude)),
                         zeros(3,length(t_vec_orbital)-1),
                         zeros(length(t_vec_orbital)),
                         zeros(3,length(t_vec_orbital)-1),
                         zeros(6,length(t_vec_attitude))
                         )


    # orbital_state = zeros(6,length(t_vec_orbital))
    # attitude_state = zeros(7,length(t_vec_attitude))
    # B_eci = zeros(3,length(t_vec_orbital)-1)
    # eclipse_hist = zeros(length(t_vec_orbital))

    B_body_sense_vec = zeros(3,length(t_vec_attitude))
    s_body_sense_vec = zeros(3,length(t_vec_attitude))
    attitude_state_sense = zeros(7,length(t_vec_attitude))

    # pre-allocate sensor stuff


    # @infiltrate
    # error()
    # initial conditions
    truth.orbital_state[:,1]  = initial_conditions.eci_rv_0
    truth.attitude_state[:,1] = [initial_conditions.ᴺqᴮ0;
                                 initial_conditions.ω0]

    # main loop
    t1 = time()
    @showprogress "Simulating..." for kk = 1:(length(t_vec_orbital)-1)

        # sun position and eclipse check
        truth.r_sun_eci[:,kk] = SD.sun_position(epc_orbital)
        eclipse = eclipse_check(truth.orbital_state[1:3,kk], truth.r_sun_eci[:,kk])
        truth.eclipse_hist[kk] = eclipse

        # atmospheric drag
        # ρ = density_harris_priester(orbital_state[1:3,kk], r_sun_eci)
        # ecef_Q_eci = SD.rECItoECEF(epc_orbital)
        # a_drag = accel_drag(orbital_state[:,kk], ρ, params.sc.mass,
        #                     params.sc.area, params.sc.cd, ecef_Q_eci)
        a_drag = zeros(3)

        # mag field vector
        truth.B_eci[:,kk] = IGRF13(truth.orbital_state[1:3,kk],epc_orbital)

        # thruster acceleration
        u_thruster = zeros(3)

        #----------------------ATTITUDE LOOP------------------------------
        # attitude dynamics inner loop
        for jj = 1:length(inner_loop_t_vec)-1

            # index of current step n
            index_n = (kk-1)*(length(inner_loop_t_vec)-1) + jj

            # true state
            ᴺqᴮ_true = truth.attitude_state[1:4,index_n]
            ω_true = truth.attitude_state[5:7,index_n]
            r_eci_true = truth.orbital_state[1:3,kk]
            v_eci_true = truth.orbital_state[4:6,kk]

            # environmental stuff
            B_eci_T_true = truth.B_eci[:,kk]
            r_sun_eci_true = truth.r_sun_eci[:,kk]
            eclipse_true = eclipse

            # ---------------truth measurements-----------------
            # attitude
            ᴺQᴮ_true = dcm_from_q(ᴺqᴮ_true)
            ᴮQᴺ_true = transpose(ᴺQᴮ_true)

            # magnetic field vector in the body
            B_body_T_true = ᴮQᴺ_true*B_eci_T_true

            # sun flux
            I_vec_true = sun_flux(r_sun_eci_true,r_eci_true,ᴺQᴮ_true)

            # sun vector in body
            s_body_true = s_body_from_I(I_vec_true)

            # eclipse sense
            eclipse_sense = eclipse_true

            # disturbance torques
            τ = zeros(3)
            sc_mag_moment = zeros(3)

            # update the attitude with RK4
            truth.attitude_state[:,index_n+1] =rk4_attitude(spacecraft_eom,
            epc_orbital, truth.attitude_state[:,index_n], sc_mag_moment, B_eci_T_true, τ, dt_attitude)

        end
        # --------------------------- end attitude loop -----------------------


        # propagate orbit one step
        truth.orbital_state[:,kk+1] =rk4_orbital(FODE, epc_orbital, truth.orbital_state[:,kk],
                                      u_thruster + a_drag, dt_orbital)

        # increment the time
        epc_orbital += dt_orbital


    end

    @show time() - t1





    ## plotting stuff
    # @infiltrate

    # plot(vec(orbital_state[1,:]),vec(orbital_state[2,:]),vec(orbital_state[3,:]))
    #
    # B_eci *=1e9
    # # plot(vec(B_eci[1,:]))
    # # plot!(vec(B_eci[2,:]))
    # # plot!(vec(B_eci[3,:]))
    # #
    # # # @infiltrate
    # t_vec_attitude = t_vec_attitude ./60
    # t_vec_orbital = t_vec_orbital ./60
    # plot(t_vec_attitude,rad2deg.(vec(attitude_state[5,:])),title = "Bdot Detumble",label = "ωₓ")
    # plot!(t_vec_attitude,rad2deg.(vec(attitude_state[6,:])))
    # plot!(t_vec_attitude,rad2deg.(vec(attitude_state[7,:])),xticks = 0:10:200,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)")
    #
    #
    # omega_norm = zeros(length(t_vec_attitude))
    # for i = 1:length(t_vec_attitude)
    #     omega_norm[i] = rad2deg(norm(vec(attitude_state[5:7,i])))
    # end
    #
    # plot(t_vec_attitude,omega_norm,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)",title = "Bdot Detumble")
    #
    # plot!(t_vec_orbital,100*eclipse_hist,label = "Eclipse")
    #
    # display(plot(vec(B_eci[1,:])))
    # display(plot!(vec(B_eci[2,:])))
    # display(plot!(vec(B_eci[3,:])))

    return sim_output = (truth=truth, t_vec_orbital = t_vec_orbital,
                         t_vec_attitude = t_vec_attitude)
end

path_to_yaml = "sim/config_attitude_test.yml"
sim_output = sim_driver(path_to_yaml)
# B_eci = sim_output.B_eci

# CSV.write("FileName.csv",  DataFrame(B_eci'), header=false)
plot(sim_output.t_vec_attitude,rad2deg.(vec(sim_output.truth.attitude_state[5,:])),title = "Bdot Detumble",label = "ωₓ")
plot!(sim_output.t_vec_attitude,rad2deg.(vec(sim_output.truth.attitude_state[6,:])))
plot!(sim_output.t_vec_attitude,rad2deg.(vec(sim_output.truth.attitude_state[7,:])),xticks = 0:10:200,xlabel = "Time (Seconds)",ylabel = "Angular Velocity (deg/s)")
