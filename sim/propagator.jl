using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator

const SD = SatelliteDynamics
const ST = SatelliteToolbox

include(joinpath(dirname(@__DIR__),"common/types.jl"))
include(joinpath(dirname(@__DIR__),"common/basis_conversions.jl"))
include(joinpath(dirname(@__DIR__),"common/misc_math_functions.jl"))
include(joinpath(dirname(@__DIR__),"dynamics/dynamics_functions.jl"))


function driver()

# load in config
include(joinpath(dirname(@__DIR__),"sim/config.jl"))


a = 7.082921757336547e6
e = 0.00069140
i = deg2rad(98.60090000)
Ω = deg2rad(127.6424)
ω= deg2rad(92.0098)
M = deg2rad(268.189)

# Declare initial state in terms of osculating orbital elements
oe0 = [a,e,i,Ω,ω,M]

# Convert osculating elements to Cartesian states
eci0 = sOSCtoCART(oe0, use_degrees=false)

# r_eci0 = eci0[1:3]
# v_eci0 = eci0[4:6]

# timing stuff
epc_orbital = SD.Epoch("2018-12-01 16:22:19.0 GPS")
dt_orbital = 10.0
tf = (24*3600)*10
t_vec = 0:dt_orbital:tf

# pre-allocate
ECI_hist = zeros(6,length(t_vec))
B_eci_hist = zeros(3,length(t_vec))

# initial conditions
ECI_hist[:,1] = eci0

# main loop
t1 = time()
for k = 1:(length(t_vec)-1)


    B_eci_hist[:,k] =


    u = zeros(3)
    ECI_hist[:,k+1] =rk4(FODE, epc_orbital, ECI_hist[:,k], u, dt_orbital)
    epc_orbital += dt_orbital

end
@show time() - t1
# @infiltrate

plot(vec(ECI_hist[1,:]),vec(ECI_hist[2,:]),vec(ECI_hist[3,:]))

end


driver()
