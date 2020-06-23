


@show @__FILE__

ss_sim_path = dirname(@__FILE__)

include(joinpath(ss_sim_path,"common/types.jl"))
include(joinpath(ss_sim_path,"common/basis_conversions.jl"))
include(joinpath(ss_sim_path,"common/misc_math_functions.jl"))
include(joinpath(ss_sim_path,"dynamics_functions.jl"))
include(joinpath(ss_sim_path,"mag_field.jl"))
include(joinpath(ss_sim_path,"bdot.jl"))
include(joinpath(ss_sim_path,"sim/config.jl"))
#
# # @infiltrate
# # error()
#
# # load in python functions
include(joinpath(ss_sim_path,"python_files/load_python_files.jl"))
