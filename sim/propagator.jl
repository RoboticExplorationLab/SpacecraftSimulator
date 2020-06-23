using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
const SD = SatelliteDynamics
const ST = SatelliteToolbox

# load in julia functions
include(joinpath(dirname(@__DIR__),"common/types.jl"))
include(joinpath(dirname(@__DIR__),"common/basis_conversions.jl"))
include(joinpath(dirname(@__DIR__),"common/misc_math_functions.jl"))
include(joinpath(dirname(@__DIR__),"dynamics/dynamics_functions.jl"))
include(joinpath(dirname(@__DIR__),"environment/mag_field.jl"))
include(joinpath(dirname(@__DIR__),"bdot.jl"))
include(joinpath(dirname(@__DIR__),"sim/config.jl"))

# load in python functions
include(joinpath(dirname(@__DIR__),"python_files/load_python_files.jl"))

function sim_driver(path_to_yaml)

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

# @infiltrate
# error()
# initial conditions
orbital_state[:,1] = initial_conditions.eci_rv_0
attitude_state[:,1] = [initial_conditions.ᴺqᴮ0;
                       initial_conditions.ω0]

# main loop
t1 = time()
for kk = 1:(length(t_vec_orbital)-1)

    # eclipse check
    r_sun_eci = SD.sun_position(epc_orbital)
    eclipse = eclipse_check(orbital_state[1:3,kk], r_sun_eci)

    @show eclipse
    # atmospheric drag
    ρ = density_harris_priester(orbital_state[1:3,kk], r_sun_eci)
    ecef_Q_eci = SD.rECItoECEF(epc_orbital)
    a_drag = accel_drag(orbital_state[:,kk], ρ, params.sc.mass, params.sc.area, params.sc.cd, ecef_Q_eci)

    # mag field vector
    B_eci[:,kk] = IGRF13(orbital_state[1:3,kk],epc_orbital)

    # thruster acceleration
    u_thruster = zeros(3)


    # attitude dynamics inner loop
    for jj = 1:length(inner_loop_t_vec)-1

        # index of current step n
        index_n = (kk-1)*(length(inner_loop_t_vec)-1) + jj

        # magnetic field vector in ECI (nT)
        B_eci_T = B_eci[:,kk]

        # bdot Control law (julia)
        m = bdot_control_law(attitude_state[1:4,index_n],
                             attitude_state[5:7,index_n],
                             params.sc.max_dipoles,
                             B_eci_T,
                             eclipse)

        # bdot control law (python)
        # m = bdot_control_law_python(attitude_state[1:4,index_n],
        #                             attitude_state[5:7,index_n],
        #                             params.sc.max_dipoles,
        #                             B_eci_T,
        #                             eclipse)

        # disturbance torques
        τ = zeros(3)

        # update the attitude with RK4
        attitude_state[:,index_n+1] =rk4_attitude(spacecraft_eom, epc_orbital, attitude_state[:,index_n], m, B_eci_T, τ, dt_attitude)

    end


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
plot(t_vec_attitude,rad2deg.(vec(attitude_state[5,:])),title = "Bdot Detumble",label = "ωₓ")
plot!(t_vec_attitude,rad2deg.(vec(attitude_state[6,:])))
plot!(t_vec_attitude,rad2deg.(vec(attitude_state[7,:])),xticks = 0:10:200,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)")


omega_norm = zeros(length(t_vec_attitude))
for i = 1:length(t_vec_attitude)
    omega_norm[i] = rad2deg(norm(vec(attitude_state[5:7,i])))
end

plot(t_vec_attitude,omega_norm,xlabel = "Time (minutes)",ylabel = "Angular Velocity (deg/s)",title = "Bdot Detumble")

end

path_to_yaml = "sim/config.yml"
sim_driver(path_to_yaml)
