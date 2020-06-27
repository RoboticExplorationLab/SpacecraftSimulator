using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
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

    # load in config from the given path
    params,initial_conditions, time_params = config(path_to_yaml)
    global params

    # timing stuff
    epc_orbital = initial_conditions.epc_orbital
    dt_orbital = time_params.dt_orbital
    dt_attitude = time_params.dt_attitude
    dt_controller = time_params.dt_controller

    # this takes care of the sample ratio
    controller_sample_ratio = Int(dt_controller/dt_attitude)

    # time vector stuff
    tf = time_params.tf
    t_vec_orbital = time_params.t_vec_orbital
    inner_loop_t_vec = time_params.inner_loop_t_vec
    t_vec_attitude = time_params.t_vec_attitude

    # pre-allocate arrays for storage
    orbital_state = zeros(6,length(t_vec_orbital))
    attitude_state = zeros(7,length(t_vec_attitude))
    B_eci = zeros(3,length(t_vec_orbital)-1)
    eclipse_hist = zeros(length(t_vec_orbital))

    B_body_sense_vec = zeros(3,length(t_vec_attitude))
    s_body_sense_vec = zeros(3,length(t_vec_attitude))
    attitude_state_sense = zeros(7,length(t_vec_attitude))

    # pre-allocate sensor stuff


    # @infiltrate
    # error()
    # initial conditions
    orbital_state[:,1] = initial_conditions.eci_rv_0
    attitude_state[:,1] = [initial_conditions.ᴺqᴮ0;
                           initial_conditions.ω0]

    # main loop
    t1 = time()
    for kk = 1:(length(t_vec_orbital)-1)

        # sun position and eclipse check
        r_sun_eci = SD.sun_position(epc_orbital)
        eclipse = eclipse_check(orbital_state[1:3,kk], r_sun_eci)
        eclipse_hist[kk] = eclipse
        # atmospheric drag
        ρ = density_harris_priester(orbital_state[1:3,kk], r_sun_eci)
        ecef_Q_eci = SD.rECItoECEF(epc_orbital)
        a_drag = accel_drag(orbital_state[:,kk], ρ, params.sc.mass,
                            params.sc.area, params.sc.cd, ecef_Q_eci)

        # mag field vector
        B_eci[:,kk] = IGRF13(orbital_state[1:3,kk],epc_orbital)

        # thruster acceleration
        u_thruster = zeros(3)


        # attitude dynamics inner loop
        attitude_loop!(orbital_state,attitude_state,kk,B_eci,r_sun_eci,
                                eclipse,inner_loop_t_vec,attitude_state_sense,
                                B_body_sense_vec,s_body_sense_vec,epc_orbital,
                                dt_attitude)

        # propagate orbit one step
        orbital_state[:,kk+1] =rk4_orbital(FODE, epc_orbital, orbital_state[:,kk],
                                      u_thruster + a_drag, dt_orbital)

        # increment the time
        epc_orbital += dt_orbital

    end

    @show time() - t1





    ## plotting stuff
    # @infiltrate

    # plot(vec(orbital_state[1,:]),vec(orbital_state[2,:]),vec(orbital_state[3,:]))
    #
    B_eci *=1e9
    # plot(vec(B_eci[1,:]))
    # plot!(vec(B_eci[2,:]))
    # plot!(vec(B_eci[3,:]))
    #
    # # @infiltrate
    t_vec_attitude = t_vec_attitude ./60
    t_vec_orbital = t_vec_orbital ./60
    plot(t_vec_attitude,rad2deg.(vec(attitude_state[5,:])),title = "Bdot Detumble",label = "ωₓ")
    plot!(t_vec_attitude,rad2deg.(vec(attitude_state[6,:])))
    plot!(t_vec_attitude,rad2deg.(vec(attitude_state[7,:])),xticks = 0:10:200,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)")


    omega_norm = zeros(length(t_vec_attitude))
    for i = 1:length(t_vec_attitude)
        omega_norm[i] = rad2deg(norm(vec(attitude_state[5:7,i])))
    end

    plot(t_vec_attitude,omega_norm,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)",title = "Bdot Detumble")

    plot!(t_vec_orbital,100*eclipse_hist,label = "Eclipse")
end

path_to_yaml = "sim/config.yml"
sim_driver(path_to_yaml)
