

struct MEKF_struct
    mu              :: Array{Array{Float64,1},1}
    Sigma           :: Array{Array{Float64,2},1}
    Q               :: Array{Float64,2} # static
    R               :: Array{Float64,2} # static
    ᴺqᴮ             :: Array{Array{Float64,1},1} # derived
    ᴺQᴮ             :: Array{Array{Float64,2},1} # derived
    ω               :: Array{Array{Float64,1},1} # derived
    β               :: Array{Array{Float64,1},1} # derived
end


function initialize_mekf_struct(time_params::NamedTuple)::MEKF_struct

    t_vec_attitude = time_params.t_vec_attitude

    attitude3 = fill(zeros(3),length(t_vec_attitude))
    attitude3x3 = fill(zeros(3,3),length(t_vec_attitude))
    attitude4 = fill(zeros(4),length(t_vec_attitude))
    attitude10 = fill(zeros(10),length(t_vec_attitude))
    attitude9x9 = fill(zeros(9,9),length(t_vec_attitude))

    MEKF = MEKF_struct(copy(attitude10),    # mu
                       copy(attitude9x9),  # Sigma
                       copy(zeros(9,9)),  # Q
                       copy(zeros(9,9)),  # R
                       copy(attitude4),    # ᴺqᴮ
                       copy(attitude3x3),  # ᴺQᴮ
                       copy(attitude3),    # ω
                       copy(attitude3))    # β


end

function initialize_mekf!(MEKF::MEKF_struct,truth::truth_state_struct)

    # process noise covariance
    MEKF.Q .= diagm(1e-12*ones(9))
    MEKF.Q[7:9,7:9] .= diagm(params.sensors.gyro.bias_noise_std^2*ones(3))

    # sensor noise covariance
    R_gyro = deg2rad(params.sensors.gyro.noise_std_degps)^2*ones(3)
    R_sun_sensor = deg2rad(params.sensors.sun_sensor.noise_std_deg)^2*ones(3)
    R_magnetometer = deg2rad(params.sensors.magnetometer.noise_std_deg)^2*ones(3)

    MEKF.R .= 1*diagm([R_gyro;R_sun_sensor;R_magnetometer])

    # initialize mu and Sigma
    MEKF.mu[1] = [truth.ᴺqᴮ[1];truth.ω[1];zeros(3)]
    MEKF.Sigma[1] = .1*eye(9)
end

function update_mekf_struct!(MEKF,index_n)

    # fill in the state
    MEKF.ᴺqᴮ[index_n] = copy(MEKF.mu[index_n][1:4])
    MEKF.ᴺQᴮ[index_n] = dcm_from_q(MEKF.ᴺqᴮ[index_n])
    MEKF.ω[index_n]   = copy(MEKF.mu[index_n][5:7])
    MEKF.β[index_n]   = copy(MEKF.mu[index_n][8:10])

end
