using LinearAlgebra, SatelliteDynamics, SatelliteToolbox

const ST = SatelliteToolbox
const SD = SatelliteDynamics


function IGRF13(r_eci,epc)

    # get decimal date
    date = SD.caldate(epc)
    year = date[1]
    day = SD.day_of_year(epc)
    decimal_date = year + day/365.2425

    # eci and ecef location
    ecef_Q_eci = SD.rECItoECEF(epc)
    eci_Q_ecef = transpose(ecef_Q_eci)

    # get ecef location
    r_ecef = ecef_Q_eci*r_eci

    # long lat geod
    longitude,latitude,altitude = SD.sECEFtoGEOC(r_ecef,use_degrees=true)

    # IGRF
    # SatelliteToolbox v0.7.1
    # B_ned_nT = igrf(decimal_date, norm(r_ecef), latitude, longitude, Val(:geocentric))
    # SatelliteToolbox v0.6.3
    # B_ned_nT = igrf12(decimal_date, norm(r_ecef), latitude, longitude, Val{:geocentric})
    B_ned_nT = my_igrf_13(decimal_date,norm(r_ecef),latitude,longitude,13)
    
    # NED and ECEF DCM
    ecef_Q_ned = ecef_Q_ned_mat(longitude,latitude)

    # conver to eci
    B_eci_nT = eci_Q_ecef*ecef_Q_ned*B_ned_nT

    # convert from nT to T
    return B_eci_nT*1e-9
end


function eclipse_check(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    """
    Computes the illumination fraction of a satellite in Earth orbit using a
    cylindrical Earth shadow model.

    Arguments:
    - `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial
                             reference frame [m; m/s]
    - `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

    Return:
    - `nu::Float64`: Illumination fraction (0 <= nu <= 1). nu = 0 means
                     spacecraft in complete shadow, nu = 1 mean spacecraft
                     fully illuminated by sun.
    References:
    1. O. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods
                                    and Applications_, 2012, p.80-83.
    """
    # Satellite position ECI
    r = x[1:3]

    # Sun-direction unit-vector
    e_sun = r_sun / norm(r_sun)

    # Projection of spacecraft position
    s = dot(r, e_sun)

    # Compute illumination
    nu = true
    if s/norm(s) >= 1.0 || norm(r - s*e_sun) > R_EARTH
        nu = false
    end

    return nu
end

function environmental_torques(truth,orb_ind,att_ind)
    """Environmental torque modeling.

    Args:
        truth: truth_state_struct
        kk: orbital index
        jj: attitude index

    Output: torque (N⋅m)

    """
    r = norm(truth.r_eci[orb_ind])
    n = normalize(transpose(truth.ᴺQᴮ[att_ind])*truth.r_eci[orb_ind])

    τ_gg = (3*GM_EARTH/(r^3))*cross(n,params.sc.J*n)

    return τ_gg
end

# function ecef_Q_ned_mat(longitude,latitude)
#
#     # ϕ = geoc[2]
#     # λ = geoc[1]
#     ϕ = latitude
#     λ = longitude
#
#     ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
#                   -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
#                    cos(ϕ)          0       -sin(ϕ)]
#
#     return ecef_Q_ned
# end
