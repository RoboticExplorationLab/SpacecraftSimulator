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
    longitude,latitude,altitude = SD.sECEFtoGEOC(r_ecef,use_degrees=false)

    # IGRF
    B_ned_nT = igrf(decimal_date, norm(r_ecef), latitude, longitude, Val(:geocentric))

    # NED and ECEF stuff
    ecef_Q_ned = ecef_Q_ned_mat(longitude,latitude)

    # conver to eci
    B_eci_nT = eci_Q_ecef*ecef_Q_ned*B_ned_nT

    return B_eci_nT
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
