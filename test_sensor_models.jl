using LinearAlgebra, SatelliteDynamics, SatelliteToolbox, Plots

using Infiltrator
using StaticArrays
const SD = SatelliteDynamics
const ST = SatelliteToolbox


# load in the julia and python functions
ss_sim_path =  dirname(@__FILE__)
include(joinpath(ss_sim_path,"load_julia_functions.jl"))

path_to_yaml = "sim/config.yml"

params,initial_conditions, time_params = config(path_to_yaml)
global params

faces = params.sc.faces

r_sun_eci = normalize([1;2;3])
r_eci = normalize([3;4;-3])




ᴺQᴮ = exp(skew_from_vec([1;2;3]))

s_body = transpose(ᴺQᴮ)*normalize(r_sun_eci - r_eci)


I_vec= sun_flux(r_sun_eci,r_eci,ᴺQᴮ)

s_body2 = s_body_from_I(I_vec)
