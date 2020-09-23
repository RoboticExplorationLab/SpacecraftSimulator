using LinearAlgebra, SatelliteDynamics



function phi_from_dcm(Q)
    # TODO: test this
    """Matrix logarithm for 3x3 orthogonal matrices (like DCM's).

    Summary:
        This is both faster and more robust than log.jl (180 degree rotations)

    Args:
        Q: orthogonal matrix (like a DCM) :: AbstractArray{Float64,2}

    Returns:
        skew symmetric matrix :: AbstractArray{Float64,2}
    """

    val = (tr(Q) - 1) / 2

    if abs(val - 1) < 1e-10
        # no rotation
        phi = [0.0; 0.0; 0.0]
    elseif abs(val + 1) < 1e-10
        # 180 degree rotation
        M = I + Q
        r = M[1, :] / norm(M[1, :])
        theta = pi
        phi = r * theta
    else
        # normal rotation (0-180 degrees)
        theta = acos(val)
        r = -(1 / (2 * sin(theta))) *
            [Q[2, 3] - Q[3, 2]; Q[3, 1] - Q[1, 3]; Q[1, 2] - Q[2, 1]]
        phi = r * theta
    end

    return phi
end


function rotz(θ)
    """Rotation matrix for rotation about the z axis"""
    return [ cos(θ)  sin(θ) 0
            -sin(θ)  cos(θ) 0
             0       0      1]
end




function myECItoECEF(seconds_from_9_1_2020_to_epoch,seconds_from_epoch)
    """ECI to ECEF rotation matrix.

    Args:
        seconds_from_9_1_2020_to_epoch: (s)
        seconds_from_epoch: (s) will be time.monotonic()

    Returns:
        ECEF_Q_ECI: DCM
    """
    t = seconds_from_9_1_2020_to_epoch + seconds_from_epoch
    return rotz(5.940294244015959 + t*7.292115146706979e-5)

    # t = seconds_from_9_1_2020_to_epoch + seconds_from_epoch
    # ang_64 = 5.940294244015959 + t*7.292115146706979e-5

    # t = convert(Float32,seconds_from_9_1_2020_to_epoch) + convert(Float32,seconds_from_epoch)
    # ang_32 = convert(Float32,5.940294244015959) + t*convert(Float32,7.292115146706979e-5)

end


epoch_0 = Epoch("2020-09-01")

epoch_anchor = Epoch("2020-09-14") -12323.3423
# test a few of these
one_month = 24*3600*30.32434 # seconds


# epc2 = epoch_0 +dt
seconds_from_9_1_2020_to_epoch = epoch_anchor - epoch_0

seconds_from_epoch = rand(-100:.00001:50)*one_month

Q = rECItoECEF(epoch_anchor + seconds_from_epoch)
Qest = myECItoECEF(seconds_from_9_1_2020_to_epoch,seconds_from_epoch )

error_phi_deg = rad2deg(norm(phi_from_dcm(Q'*Qest)))

@show error_phi_deg
@show max_error_distance_km = 1.08*R_EARTH*2*pi*(error_phi_deg/360)/1000

@show Q'*Qest
# theta_0 = 5.940294244015959
#
# omega = 7.292115146706979e-5
#
# Q_true = rECItoECEF(epc)
#
# Q_est = rotz(theta_0)
#
# dt = -1023442
# new_theta = theta_0 + omega*dt
#
# Q_est2 = rotz(new_theta)
#
# Q_true2 = rECItoECEF(epc + dt)


# now lets try some 6 digit stuff

# 7.292115146706979e-5

t = 5.234*one_month
w_earth_pt1 = convert(Float32,7.29211e-5)
w_earth_pt2 = convert(Float32,.000005146706e-5)

t_32 = 5.234 * convert(Float32,one_month)

ang_64 = t*OMEGA_EARTH % 2*pi

ang_32 = t_32*w_earth_pt1 % (2*convert(Float32,pi))
ang_32 = t_32*w_earth_pt1 % (2*convert(Float32,pi)) + t_32*w_earth_pt2 % (2*convert(Float32,pi))
ang_32 = t_32*w_earth_pt1 % 2*pi


# julia is awesome
ᴺRᶜ = ᴺRᴬ * ᴬRᴮ * ᴮRᶜ

# this is still ok
R_n_c = R_n_a * R_a_b * R_b_c

# atrocious
R_c_to_n = R_a_to_n * R_b_to_a * R_c_to_b


end
