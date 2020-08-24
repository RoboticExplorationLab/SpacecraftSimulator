using LinearAlgebra, SatelliteToolbox

R_EARTH     = 6.378136300e6 # m
decimal_date = 2020.8
λ  = deg2rad(35)
Ω = deg2rad(15)
r_norm = 1.05*R_EARTH


B1 = igrf(decimal_date, r_norm, λ, Ω, Val(:geocentric))
B2 = igrf13syn(0,decimal_date,2,r_norm/1000,colatd_from_latd(rad2deg(λ)), rad2deg(Ω))
B3 = my_igrf(gh_igrf13_trim,decimal_date,r_norm/1000,rad2deg(λ),rad2deg(Ω))
# B4 = my_igrf(gh_igrf13,decimal_date,r_norm/1000,colatd_from_latd(rad2deg(λ)), rad2deg(Ω))
B1 = [B1[1];B1[2];B1[3]]
B2 = [B2[1];B2[2];B2[3]]


@show B1-B2
@show B2-B3
# B1 = igrf12(decimal_date, r_norm, λ, Ω, Val{:geocentric})
#
# B2 = igrf12syn(0,decimal_date,2,r_norm, 90 - rad2deg(λ), rad2deg(Ω)+180)

# B2 = [B2[1];B2[2];B2[3]]*1e9
#
# @show B1
# @show B2

# colatitude should be
# B1 = igrf(decimal_date,[])

function colatd_from_latd(latd)

    return 90 -latd
end
