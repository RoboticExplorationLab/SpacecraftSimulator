
ss_sim_path =  dirname(dirname(@__FILE__))

cd(joinpath(ss_sim_path,"virtual_env"))
Pkg.activate(".")


function run_tests()
    include(joinpath(ss_sim_path,"test/attitude_fx_tests.jl"))
    include(joinpath(ss_sim_path,"test/sensor_tests.jl"))
    include(joinpath(ss_sim_path,"test/igrf_tests.jl"))
    include(joinpath(ss_sim_path,"test/propagator_reg_tests.jl"))

end
