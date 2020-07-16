

#------------------------------TRUTH models------------------------------------

function sun_body_normalized(r_sun_eci,r_eci,ᴺQᴮ)
    # posiiton vector from spacecraft to sun
    sc_r_sun = r_sun_eci - r_eci

    # normalize and express in the body frame
    return transpose(ᴺQᴮ)*normalize(sc_r_sun)
end

function sun_flux(r_sun_eci,r_eci,ᴺQᴮ,eclipse)
    """Gets the true flux values for each coarse sun sensor

    Args:
        r_sun_eci: position vector from eci to sun
        r_eci: sc position vector from eci to sc
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

function s_body_from_I(I_vec)
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




function create_LEM_offset(scale_factor,cross_axis_sensitivity)

    scale_mat = sqrt(diagm(scale_factor))*randn(length(scale_factor))

    cas = cross_axis_sensitivity
    cross_axis_mat = [1                  sqrt(cas)*randn() sqrt(cas)*randn();
                      sqrt(cas)*randn()  1                 sqrt(cas)*randn();
                      sqrt(cas)*randn()  sqrt(cas)*randn() 1                ]

    A_error = (I + scale_mat)*cross_axis_mat

    return A_error

end

function add_LEM_noise(true_vector,A_error, bias, sqrt_noise_cov)

    return A_error*true_vector + bias + sqrt_noise_cov*randn(length(true_vector))

end
