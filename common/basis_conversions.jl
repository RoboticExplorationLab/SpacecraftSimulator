using LinearAlgebra, SatelliteDynamics

const SD = SatelliteDynamics

function ecef_Q_ned_mat(longitude::Float64,latitude::Float64)::Matrix

    # ϕ = geoc[2]
    # λ = geoc[1]
    ϕ = latitude
    λ = longitude

    ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
                  -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
                   cos(ϕ)          0.0     -sin(ϕ)]

    return ecef_Q_ned
end

function geoc_from_ecef(r_ecef::Vector)::Tuple{Float64,Float64}
    x,y,z = ecef

    lat = atan(z, sqrt(x*x + y*y))
    lon = atan(y, x)

    return lat, lon
end

# lon = pi/3
# lat = -pi/4
#
# ecef_Q_ned = ecef_Q_ned_mat(lon,lat)
#
# # Compute ENZ basis vectors
#     eE = [-sin(lon) ; cos(lon) ; 0]
#     eN = [-sin(lat)*cos(lon) ; -sin(lat)*sin(lon) ; cos(lat)]
#     eZ = [cos(lat)*cos(lon) ; cos(lat)*sin(lon) ; sin(lat)]
#
#     # Construct Rotation matrix
#     enz_Q_ecef = hcat(eE, eN, eZ)'
#
# ecef_Q_enz = enz_Q_ecef'
#
#
# N1 = ecef_Q_ned[:,1]
# E1 = ecef_Q_ned[:,2]
# D1 = ecef_Q_ned[:,3]
#
# N2 = ecef_Q_enz[:,2]
# E2 = ecef_Q_enz[:,1]
# D2 = -ecef_Q_enz[:,3]
