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

    # conver to eci
    B_eci_nT = eci_Q_ecef*ecef_Q_ned*B_ned_nT

    return B_eci_nT
end




function ecef_Q_ned_mat(longitude,latitude)

    # ϕ = geoc[2]
    # λ = geoc[1]
    ϕ = latitude
    λ = longitude

    ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
                  -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
                   cos(ϕ)          0       -sin(ϕ)]

    return ecef_Q_ned
end

# start the epoch
epc = SD.Epoch("2018-12-01 16:22:19.0 GPS")

# random eci position vector
r_eci = 1.1*R_EARTH*normalize(randn(3))

# eci and ecef transforms
ecef_Q_eci = SD.rECItoECEF(epc)
eci_Q_ecef = transpose(ecef_Q_eci)

# get ecef location
r_ecef = ecef_Q_eci*r_eci

# geocentric coordinates
# geoc = SD.sECEFtoGEOC(r_ecef,use_degrees=false)
longitude,latitude,altitude = SD.sECEFtoGEOC(r_ecef,use_degrees=false)
# longitude
# latitude
# altitude


# longitude = geoc[1]
# latitude = geoc[2]
# altitude = geoc[3]
#
# ϕ = geoc[2]
# λ = geoc[1]
#
# ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
#               -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
#                cos(ϕ)          0       -sin(ϕ)]
#
# n = [1;0;0]
# e = [0;1;0]
# d = [0;0;1]
#
# ecef_Q_ned*d


B_ned_nT = igrf(2017.12313, norm(r_ecef), latitude, longitude, Val(:geocentric))


# B_eci_nT = eci_Q_ecef*ecef_Q_ned*B_ned_nT


function ecef_Q_ned_mat(longitude,latitude)

    # ϕ = geoc[2]
    # λ = geoc[1]
    ϕ = latitude
    λ = longitude

    ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
                  -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
                   cos(ϕ)          0       -sin(ϕ)]

    return ecef_Q_ned
end






r_sun = sun_position(epc)

r_sun = norm(r_sun)*[1;0;0]
# Satellite inertial position
r = r_eci
r = norm(r_eci)*[0;1;0]

# Sun-direction unit-vector
e_sun = r_sun / norm(r_sun)

# Projection of spacecraft position
s = dot(r, e_sun)

# Compute illumination
nu = 0.0
if s/norm(s) >= 1.0 || norm(r - s*e_sun) > R_EARTH
    nu = 1.0
end

return nu






# igrf(2017.12313, alt, lat, long, Val(:geodetic))
# use_degrees = false
#
# # duncan's part
# # Extract lat and lon
# lon = geoc[1]
# lat = geoc[2]
#
# # Handle non-explict use-degrees
# if length(geoc) == 3
#     alt = geoc[3]
# else
#     alt = 0.0
# end
#
# # Convert input to radians
# if use_degrees
#     lat = lat*pi/180.0
#     lon = lon*pi/180.0
# end
#
# # Check validity of input
# if lat < -pi/2 || lat > pi/2
#     throw(ArgumentError("Lattiude, $lat, out of range. Must be between -90 and 90 degrees."))
# end
#
# # Compute Earth fixed coordinates
# r       = WGS84_a + alt
# x = r*cos(lat)*cos(lon)
# y = r*cos(lat)*sin(lon)
# z = r*sin(lat)
