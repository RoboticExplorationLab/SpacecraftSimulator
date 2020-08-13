# @show @__FILE__

ss_sim_path = dirname(@__FILE__)

include(joinpath(ss_sim_path,"common/types.jl"))
include(joinpath(ss_sim_path,"common/basis_conversions.jl"))
include(joinpath(ss_sim_path,"common/misc_math_functions.jl"))
include(joinpath(ss_sim_path,"dynamics_functions.jl"))
include(joinpath(ss_sim_path,"mag_field.jl"))
include(joinpath(ss_sim_path,"bdot.jl"))
include(joinpath(ss_sim_path,"sim/config.jl"))
include(joinpath(ss_sim_path,"common/testing_functions.jl"))
include(joinpath(ss_sim_path,"common/state_structs.jl"))
include(joinpath(ss_sim_path,"sensor_models.jl"))
include(joinpath(ss_sim_path,"inner_attitude_loop.jl"))
include(joinpath(ss_sim_path,"MEKF/MEKF_functions.jl"))
include(joinpath(ss_sim_path,"MEKF/MEKF_utils.jl"))
#
# # @infiltrate
# # error()
#
# # load in python functions TODO: fix this with brian tracy
# include(joinpath(ss_sim_path,"python_files/load_python_files.jl"))
