

#------------------------------TRUTH models------------------------------------

function sun_body_normalized(r_sun_eci::Vec,r_eci::Vec,ᴺQᴮ::Mat)::Vec
    # positon vector from spacecraft to sun
    sc_r_sun = r_sun_eci - r_eci

    # normalize and express in the body frame
    return transpose(ᴺQᴮ)*normalize(sc_r_sun)
end

function sun_flux(r_sun_eci::Vec,r_eci::Vec,ᴺQᴮ::Mat,eclipse::Bool)::Vec
    """Gets the true flux values for each coarse sun sensor

    Args:
        r_sun_eci: position vector from eci0 to sun
        r_eci: sc position vector from eci0 to sc
        ᴺQᴮ: DCM relating eci (N) and body (B)
        eclipse:

    Returns:
        I_vec: flux values on each face (0-1)
    """

    if eclipse
        return zeros(6)
    else
        # normalize and express in the body frame
        s = sun_body_normalized(r_sun_eci,r_eci,ᴺQᴮ)

        # Get dot products
        I_vec = params.sc.faces*s

        # only keep the positive ones (negative ones go to 0)
        I_vec = I_vec .* (I_vec .> 0.0)

        return I_vec
    end

end

function s_body_from_I(I_vec::Vec)::Vec
    """Get unit sun vector expressed in the body frame from solar flux values.

    Args:
        I_vec: flux values on each face (0-1)

    Returns:
        s_body: unit vector from sc to sun expressed in body
    """

    s_body = normalize([(I_vec[1]-I_vec[2]);
                        (I_vec[3]-I_vec[4]);
                        (I_vec[5]-I_vec[6])])

    return s_body
end


#------------------Linear Error Model------------------------------------------


function sample_inertia(J::Mat,deg_std::Float64,scale_std::Float64)::Mat

    # take eigen decomposition
    eigen_decomp = eigen(J)
    D = diagm(eigen_decomp.values)
    S = eigen_decomp.vectors

    # scale the moments
    D_sample = D*diagm(ones(3)+ scale_std*randn(3))

    # rotate them
    R_sample = exp(hat(deg2rad(deg_std)*randn(3)))
    J_sample = R_sample*S*D_sample*S'*R_sample'

    return J_sample
end

function S03_noise(vector::Vec, noise_std::Float64)
    # here we add noise S03 style

    # noise axis angle vector
    phi_noise = noise_std*randn(3)

    # apply the noise rotation to the vector
    noisy_vector = dcm_from_phi(phi_noise)*vector

    return noisy_vector
end

function measurements(truth::truth_state_struct,orb_ind::Int,index_n::Int)

    # genereate gyro
    true_ω = truth.ω[index_n]
    gyro_offset = params.sensors.offsets.gyro
    # TODO: include bias dynamics in there somewhere
    # gyro_bias = params.sensors.offsets.gyro_bias
    gyro_bias = truth.β[index_n]
    gyro_noise = deg2rad(params.sensors.gyro.noise_std_degps)*randn(3)
    meas_ω = gyro_offset*true_ω + gyro_bias + gyro_noise

    # generate sun sensor
    true_sun_body = truth.sun_body[index_n]
    sun_sensor_offset = params.sensors.offsets.sun_sensor
    meas_sun_sensor = S03_noise(sun_sensor_offset*true_sun_body,
                      deg2rad(params.sensors.sun_sensor.noise_std_deg))

    # generate magnetometer
    true_B_body = normalize(truth.B_body[index_n])
    magnetometer_offset = params.sensors.offsets.magnetometer
    meas_B_body = S03_noise(magnetometer_offset*true_B_body,
                      deg2rad(params.sensors.magnetometer.noise_std_deg))

    # @infiltrate
    # error()
    return meas_ω, meas_sun_sensor, meas_B_body

end






# function create_LEM_offset(scale_factor::Vec,cross_axis_sensitivity::Float64)::Mat
#
#     scale_mat = sqrt(diagm(scale_factor))*randn(length(scale_factor))
#
#     cas = cross_axis_sensitivity
#     cross_axis_mat = [1                  sqrt(cas)*randn() sqrt(cas)*randn();
#                       sqrt(cas)*randn()  1                 sqrt(cas)*randn();
#                       sqrt(cas)*randn()  sqrt(cas)*randn() 1                ]
#
#     A_error = (I + scale_mat)*cross_axis_mat
#
#     return A_error
#
# end
#
# function add_LEM_noise(true_vector,A_error, bias, sqrt_noise_cov)
#
#     return A_error*true_vector + bias + sqrt_noise_cov*randn(length(true_vector))
#
# end
