using LinearAlgebra


function rotz(θ)
    """Rotation matrix for rotation about the z axis"""
    return [ cos(θ)  sin(θ) 0
            -sin(θ)  cos(θ) 0
             0       0      1]
end


epoch_0 = Epoch("2020-09-01")

function myECItoECEF(seconds_since_9_1_2020)
    return rotz(5.940294244015959 + seconds_since_9_1_2020*7.292115146706979e-5)
end


# test a few of these
one_month = 24*3600*30 # seconds

dt = - 2*one_month
# epc2 = epoch_0 +dt

Q = rECItoECEF(epoch_0 + dt)
Qest = myECItoECEF(dt)

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
