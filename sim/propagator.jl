using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
const SD = SatelliteDynamics
const ST = SatelliteToolbox

include(joinpath(dirname(@__DIR__),"common/types.jl"))
include(joinpath(dirname(@__DIR__),"common/basis_conversions.jl"))
include(joinpath(dirname(@__DIR__),"common/misc_math_functions.jl"))
include(joinpath(dirname(@__DIR__),"dynamics/dynamics_functions.jl"))
include(joinpath(dirname(@__DIR__),"environment/mag_field.jl"))

function driver()

# load in config
include(joinpath(dirname(@__DIR__),"sim/config.jl"))


# a = 7.082921757336547e6
# e = 0.00069140
# i = deg2rad(98.60090000)
# Ω = deg2rad(127.6424)
# ω= deg2rad(92.0098)
# M = deg2rad(268.189)
a = 6.7780e6
e = 0
i = .6102
Ω = 1.4175
ω= deg2rad(92.0098)
M = 2.088

# Declare initial state in terms of osculating orbital elements
oe0 = [a,e,i,Ω,ω,M]

# Convert osculating elements to Cartesian states
eci0 = sOSCtoCART(oe0, use_degrees=false)

# r_eci0 = eci0[1:3]
# v_eci0 = eci0[4:6]

# timing stuff
epc_orbital = SD.Epoch("2018-12-01 16:22:19.0 GPS")
epc_orbital = SD.Epoch("2011-05-10 04:56:19.0 GPS")
dt_orbital = 10.0
dt_attitude = .01
dt_controller = .1

# this takes care of the sample ratio
controller_sample_ratio = Int(dt_controller/dt_attitude)

# time vector stuff
tf = (24*3600)*.1
t_vec_orbital = 0:dt_orbital:tf
inner_loop_t_vec = 0:dt_attitude:dt_orbital
t_vec_attitude = 0:dt_attitude:(t_vec_orbital[end]+dt_orbital-dt_attitude)

# pre-allocate
orbital_state = zeros(6,length(t_vec_orbital))
attitude_state = zeros(7,Int(length(t_vec_orbital)*dt_orbital*(1/dt_attitude)))
B_eci = zeros(3,length(t_vec_orbital)-1)

# initial conditions
orbital_state[:,1] = eci0
q0 = [0;0;0;1]
# ω0 = [.2;7;.2]
ω0 = deg2rad(360)*normalize(randn(3))
ω0 = [.01;.01;.01]
attitude_state[:,1] = [q0;ω0]

# main loop
t1 = time()
for kk = 1:(length(t_vec_orbital)-1)

    # eclipse check
    r_sun_eci = SD.sun_position(epc_orbital)
    eclipse = eclipse_check(orbital_state[1:3,kk], r_sun_eci)

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

        # magnetic moment (control law)
        # m = zeros(3)
        q = attitude_state[1:4,index_n]
        n_Q_b = dcm_from_q(q)
        ω = attitude_state[5:7,index_n]

        bdot = -cross(ω,n_Q_b'*B_eci_T)

        # m = -.5*params.sc.max_dipoles .* sign.(bdot)


        m = -5*bdot/norm(B_eci_T)*1e4
        # @show m

        # disturbance torques
        τ = zeros(3)


        # @infiltrate
        # error()


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
plot(t_vec_attitude,vec(attitude_state[5,:]))
plot!(t_vec_attitude,vec(attitude_state[6,:]))
plot!(t_vec_attitude,vec(attitude_state[7,:]),xticks = 0:1:100)

end


driver()
