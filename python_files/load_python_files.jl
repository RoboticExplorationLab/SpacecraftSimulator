@show pwd()

using PyCall

pushfirst!(PyVector(pyimport("sys")."path"), "")
pushfirst!(PyVector(pyimport("sys")."path"), "python_files")

a = pyimport("python_files")."test1"

b = pyimport("python_files")."test2"
