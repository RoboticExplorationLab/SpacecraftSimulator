using LinearAlgebra, SatelliteToolbox

# R_EARTH     = 6.378136300e6 # m
# decimal_date = 2020.8
# λ  = deg2rad(35)
# Ω = deg2rad(15)
# r_norm = 1.05*R_EARTH

function rand_in_range(lower::Real, upper::Real)::Real
    """Random number within range with a uniform distribution.

    Args:
        lower: lower value in range
        upper: upper value in range

    Returns:
        random number in the specified range
    """

    delta = upper - lower
    return lower + rand() * delta
end
function colatd_from_latd(latd)

    return 90 -latd
end

# test the stuff


function ig_test(order)

    N = 10000
    ang_err = zeros(N)
    mag_err = zeros(N)
for i = 1:N

    R_EARTH     = 6.378136300e6 # m
    decimal_date = 2020.8
    λ  = deg2rad(rand_in_range(-90,90))
    Ω = deg2rad(rand_in_range(-180,180))
    r_norm = 1.05*R_EARTH


B1 = igrf(decimal_date, r_norm, λ, Ω, Val(:geocentric))
B2 = igrf13syn(0,decimal_date,2,r_norm/1000,colatd_from_latd(rad2deg(λ)), rad2deg(Ω))
B3 = my_igrf(gh_igrf13_trim,decimal_date,r_norm/1000,rad2deg(λ),rad2deg(Ω),13)
# B4 = my_igrf(gh_igrf13,decimal_date,r_norm/1000,colatd_from_latd(rad2deg(λ)), rad2deg(Ω))
B1 = [B1[1];B1[2];B1[3]]
B2 = [B2[1];B2[2];B2[3]]


if norm(B1-B2)>1e-3
    error("first failed")
elseif norm(B2-B3)>1e-3
    error("second failed")
end


# now let's check the 6th order stuff
B6 = my_igrf(gh_igrf13_trim,decimal_date,r_norm/1000,rad2deg(λ),rad2deg(Ω),order)

ang_err[i] = acos(dot(normalize(B3),normalize(B6)))
nb3 = norm(B3)
nb6 = norm(B6)

mag_err[i] = 100*(nb6-nb3)/nb3
# @show norm(B2-B1)
end

return ang_err, mag_err

end

order = 4
ang_err, mag_err = ig_test(order)


mat"
figure
hold on
title(sprintf('Angle Error for %dth Order IGRF',$order))
histogram(rad2deg($ang_err))
xlabel('Angle Error (deg)')
ylabel('Frequency')
hold off
"

mat"
figure
hold on
title(sprintf('Magnitude Error for %dth Order IGRF',$order))
histogram(($mag_err))
xlabel('Magnitude Error (%)')
ylabel('Frequency')
hold off
"
# B1 = igrf12(decimal_date, r_norm, λ, Ω, Val{:geocentric})
#
# B2 = igrf12syn(0,decimal_date,2,r_norm, 90 - rad2deg(λ), rad2deg(Ω)+180)

# B2 = [B2[1];B2[2];B2[3]]*1e9
#
# @show B1
# @show B2

# colatitude should be
# B1 = igrf(decimal_date,[])
