"""Functions relating to the state struct."""

mutable struct truth_state_struct
    r_eci           :: Array{Array{Float64,1},1}
    v_eci           :: Array{Array{Float64,1},1}
    ᴺqᴮ             :: Array{Array{Float64,1},1}
    ᴺQᴮ             :: Array{Array{Float64,2},1} # derived
    ω               :: Array{Array{Float64,1},1}
    B_eci           :: Array{Array{Float64,1},1} # derived
    eclipse_hist    :: Array{Bool,1}             # derived
    r_sun_eci       :: Array{Array{Float64,1},1} # derived
    I_sun_flux      :: Array{Array{Float64,1},1} # derived
    sun_body        :: Array{Array{Float64,1},1} # derived
    B_body          :: Array{Array{Float64,1},1} # derived
    epc_orbital     :: Epoch
    β               :: Array{Array{Float64,1},1}
end

# struct MEKF_struct
#     mu              :: Array{Array{Float64,1},1}
#     Sigma           :: Array{Array{Float64,2},1}
#     Q               :: Array{Float64,2} # static
#     R               :: Array{Float64,2} # static
#     ᴺqᴮ             :: Array{Array{Float64,1},1} # derived
#     ᴺQᴮ             :: Array{Array{Float64,2},1} # derived
#     ω               :: Array{Array{Float64,1},1} # derived
#     β               :: Array{Array{Float64,1},1} # derived
# end

struct sensor_state_struct
    ω               :: Array{Array{Float64,1},1}
    I_sun_flux      :: Array{Array{Float64,1},1}
    sun_body        :: Array{Array{Float64,1},1} # derived
    B_body          :: Array{Array{Float64,1},1}
    J               :: Array{Float64,2}
    invJ            :: Array{Float64,2}
    sun_eci         :: Array{Array{Float64,1},1}
    B_eci           :: Array{Array{Float64,1},1}
end

function initialize_sensors_struct(time_params::NamedTuple)::sensor_state_struct

        """Initialize the sensor struct with pre-allocated arrays and IC's."""
        t_vec_orbital  = time_params.t_vec_orbital
        t_vec_attitude = time_params.t_vec_attitude

        orbital3 = fill(zeros(3),length(t_vec_orbital))
        attitude3 = fill(zeros(3),length(t_vec_attitude))
        attitude4 = fill(zeros(4),length(t_vec_attitude))
        attitude3x3 = fill(zeros(3,3),length(t_vec_attitude))
        attitude6 = fill(zeros(6),length(t_vec_attitude))

        sensors = sensor_state_struct(copy(attitude3),      # ω
                                      copy(attitude6),      # I_sun_flux
                                      copy(attitude3),      # sun_body
                                      copy(attitude3),      # B_body
                                      params.sensors.J,     # J
                                      params.sensors.invJ,   # invJ
                                      copy(orbital3),        # sun_eci
                                      copy(orbital3))        # B_eci
        return sensors
end

# struct MEKF_struct
#     mu              :: Array{Array{Float64,1},1}
#     Sigma           :: Array{Array{Float64,2},1}
#     Q               :: Array{Float64,2} # static
#     R               :: Array{Float64,2} # static
#     ᴺqᴮ             :: Array{Array{Float64,1},1} # derived
#     ᴺQᴮ             :: Array{Array{Float64,2},1} # derived
#     ω               :: Array{Array{Float64,1},1} # derived
#     β               :: Array{Array{Float64,1},1} # derived
# end

# function initialize_mekf_struct(time_params::NamedTuple)::MEKF_struct
#
#     t_vec_attitude = time_params.t_vec_attitude
#
#     attitude3 = fill(zeros(3),length(t_vec_attitude))
#     attitude3x3 = fill(zeros(3,3),length(t_vec_attitude))
#     attitude4 = fill(zeros(4),length(t_vec_attitude))
#     attitude10 = fill(zeros(10),length(t_vec_attitude))
#     attitude9x9 = fill(zeros(9,9),length(t_vec_attitude))
#
#     MEKF = MEKF_struct(copy(attitude10),    # mu
#                        copy(attitude9x9),  # Sigma
#                        copy(zeros(9,9)),  # Q
#                        copy(zeros(9,9)),  # R
#                        copy(attitude4),    # ᴺqᴮ
#                        copy(attitude3x3),  # ᴺQᴮ
#                        copy(attitude3),    # ω
#                        copy(attitude3))    # β
#
#
# end


function sensors_update!(sensors, truth,orb_ind,index_n)

    meas_ω, meas_sun_sensor, meas_B_body = measurements(truth,orb_ind,index_n)

    sensors.ω[index_n] = meas_ω
    sensors.sun_body[index_n] = meas_sun_sensor
    sensors.B_body[index_n] = meas_B_body

    # sun and B vectors in the body
    ᴺQᴮ = truth.ᴺQᴮ[index_n]
    sensors.sun_eci[orb_ind] = ᴺQᴮ * truth.sun_body[index_n]
    sensors.B_eci[orb_ind] = truth.B_eci[orb_ind]

    # @infiltrate
    # error()

end

# function initialize_mekf!(MEKF,truth)
#
#     MEKF.Q .= diagm(1e-8*ones(9))
#     MEKF.Q[7:9,7:9] .= diagm(1e-3*ones(3))
#
#     R_gyro = deg2rad(params.sensors.gyro.noise_std_degps)^2*ones(3)
#     R_sun_sensor = deg2rad(params.sensors.sun_sensor.noise_std_deg)^2*ones(3)
#     R_magnetometer = deg2rad(params.sensors.magnetometer.noise_std_deg)^2*ones(3)
#
#     MEKF.R .= 1.3*diagm([R_gyro;R_sun_sensor;R_magnetometer])
#
#     # initialize mu and Sigma
#     MEKF.mu[1] = [truth.ᴺqᴮ[1];truth.ω[1];zeros(3)]
#     MEKF.Sigma[1] = .1*eye(9)
# end
#
# function update_mekf_struct!(MEKF,index_n)
#
#     # fill in the state
#     MEKF.ᴺqᴮ[index_n] = copy(MEKF.mu[index_n][1:4])
#     MEKF.ᴺQᴮ[index_n] = dcm_from_q(MEKF.ᴺqᴮ[index_n])
#     MEKF.ω[index_n]   = copy(MEKF.mu[index_n][5:7])
#     MEKF.β[index_n]   = copy(MEKF.mu[index_n][8:10])
#
# end

function orbital_truth_struct_update!(truth::truth_state_struct,k::Int)
    """Fill in the truth struct with derived properties for orbital changes."""

    # sun position
    truth.r_sun_eci[k] = SD.sun_position(truth.epc_orbital)

    # eclipse
    truth.eclipse_hist[k] = eclipse_check(truth.r_eci[k],truth.r_sun_eci[k])

    # ECI magnetic field (T)
    truth.B_eci[k] = IGRF13(truth.r_eci[k],truth.epc_orbital)

end

function attitude_truth_struct_update!(truth::truth_state_struct,k::Int,jj::Int)
    """Fill in the truth struct with derived properties for attitude changes."""

    # pull attitude and get the DCM
    truth.ᴺQᴮ[jj]   = dcm_from_q(truth.ᴺqᴮ[jj])

    # magnetic field in the body
    truth.B_body[jj] = transpose(truth.ᴺQᴮ[jj])*truth.B_eci[k]

    # Sun Flux
    truth.I_sun_flux[jj] = sun_flux(truth.r_sun_eci[k],truth.r_eci[k],
                                   truth.ᴺQᴮ[jj],truth.eclipse_hist[k])

    # sun_body
    truth.sun_body[jj] = sun_body_normalized(truth.r_sun_eci[k],truth.r_eci[k],
                                            truth.ᴺQᴮ[jj])

end


function initialize_struct(struct_type_name   ::DataType,
                           time_params        ::NamedTuple,
                           initial_conditions ::NamedTuple)

    """Initialize the truth struct with pre-allocated arrays and IC's."""
    t_vec_orbital  = time_params.t_vec_orbital
    t_vec_attitude = time_params.t_vec_attitude

    orbital3 = fill(zeros(3),length(t_vec_orbital))
    orbital6 = fill(zeros(6),length(t_vec_orbital))
    attitude3 = fill(zeros(3),length(t_vec_attitude))
    attitude4 = fill(zeros(4),length(t_vec_attitude))
    attitude3x3 = fill(zeros(3,3),length(t_vec_attitude))
    attitude6 = fill(zeros(6),length(t_vec_attitude))
    attitude7 = fill(zeros(7),length(t_vec_attitude))
    orbital_bool = fill(false,length(t_vec_orbital))
    epc_orbital = Epoch("2018-12-20")

    truth = struct_type_name(copy(orbital3),      # r_eci
                             copy(orbital3),      # v_eci
                             copy(attitude4),     # q
                             copy(attitude3x3),   # Q
                             copy(attitude3),     # ω
                             copy(orbital3),      # B_eci
                             copy(orbital_bool),  # eclipse
                             copy(orbital3),      # r_sun_eci
                             copy(attitude6),     # I_sun_flux
                             copy(attitude3),     # sun_body
                             copy(attitude3),     # B_body
                             epc_orbital,         # time
                             copy(attitude3)      # β
                             )

    # initial conditions
    truth.r_eci[1]  = initial_conditions.eci_rv_0[1:3]
    truth.v_eci[1]  = initial_conditions.eci_rv_0[4:6]
    truth.ᴺqᴮ[1] =  initial_conditions.ᴺqᴮ0
    truth.ω[1] = initial_conditions.ω0
    truth.epc_orbital = initial_conditions.epc_orbital
    truth.β[1] = params.sensors.offsets.gyro_bias

        return truth
end
