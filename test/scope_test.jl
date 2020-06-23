using LinearAlgebra


function driver()

    include(joinpath(pwd(),"sim/config.jl"))

    a = 3

    b = subfx()

    @show a
    @show b

end

function subfx()

    b = 3

    @show J
end

driver()
