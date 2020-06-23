@show pwd()

using PyCall

include(joinpath(dirname(@__DIR__),"common/types.jl"))
include(joinpath(dirname(@__DIR__),"common/basis_conversions.jl"))
include(joinpath(dirname(@__DIR__),"common/misc_math_functions.jl"))
include(joinpath(dirname(@__DIR__),"dynamics/dynamics_functions.jl"))
include(joinpath(dirname(@__DIR__),"environment/mag_field.jl"))
include(joinpath(dirname(@__DIR__),"bdot.jl"))

if PyVector(pyimport("sys")."path")[2] != "python_files"
    pushfirst!(PyVector(pyimport("sys")."path"), "")
    pushfirst!(PyVector(pyimport("sys")."path"), "python_files")
end


# a = pyimport("python_files")."test1"

# b = pyimport("python_files")."test2"

# c = pyimport("python_files")."foo"

bdot_control_law_python = pyimport("python_files")."bdot_control_law"


# q = normalize(randn(4))
# w = randn(3)
# max_dipoles = rand(3)
# B_eci_T = randn(3)
# eclipse = false
#
# m = bdot_control_law_python(q,w,max_dipoles,B_eci_T,eclipse)
# m2 = bdot_control_law(q,w,max_dipoles,B_eci_T,eclipse)
# @show m
# @show m2
