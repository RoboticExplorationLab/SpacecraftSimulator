

mutable struct truth_struct
    orbital_state   :: Array{Float64,2}
    attitude_state  :: Array{Float64,2}
    ᴺqᴮ             :: Array{Float64,1}
    ᴺQᴮ             :: Array{Float64,2}
    ω               :: Array{Float64,1}
    B_eci           :: Array{Float64,2}
    eclipse_hist    :: Array{Bool,1}
    r_sun_eci       :: Array{Float64,2}
    I_sun_flux      :: Array{Float64,2}
    sun_body        :: Array{Float64,2}
    B_body          :: Array{Float64,2}
end

mutable struct sense_struct
    orbital_state   :: Array{Float64,2}
    attitude_state  :: Array{Float64,2}
    B_eci           :: Array{Float64,2}
    B_body          :: Array{Float64,2}
    eclipse_hist    :: Array{Bool,1}
end
