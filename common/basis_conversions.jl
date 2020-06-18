using LinearAlgebra, SatelliteDynamics

const SD = SatelliteDynamics

function ecef_Q_ned_mat(longitude,latitude)

    # ϕ = geoc[2]
    # λ = geoc[1]
    ϕ = latitude
    λ = longitude

    ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
                  -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
                   cos(ϕ)          0.0     -sin(ϕ)]

    return ecef_Q_ned
end
