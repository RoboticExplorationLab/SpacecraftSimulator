using LinearAlgebra, SatelliteDynamics, SOFA, Attitude

# epc = Epoch(2020,11,9,20,31,4,0,tsys = "UTC")
#
# jd_current = jd(epc)
#
# # θ = earth_rotation(epc)
#
# TU = (jd_current - 2451545)
#
# θ2 = 2*pi*0.7790572732640 + 2*pi*1.00273781191135448*TU
#
# theta = wrap_to_2pi(θ2)
#
#
# epc = Epoch(2020,11,9,20,31,4,0,tsys = "UTC") - 86164.09890369035*20
#
# # jd_current = jd(epc)
#
# # θ = earth_rotation(epc)
# era = iauEra00(MJD_ZERO, mjd(epc, tsys="UT1"))
#
# mjd1 = 59161
#
# mjd2 = 1.85491
#
#
# earth_period = 2*pi/OMEGA_EARTH

epc = Epoch(2020,11,9,2,34,32,0,tsys = "UTC")

# mjd_current = mjd(epc)
mjd_current = 59163 + 0.204329

# θ = earth_rotation(epc)

# TU = (mjd_current -51544.5)

# θ2 = 2*pi*0.7790572732640 + 2*pi*1.00273781191135448*mjd_current - 324750.3228110344
θ2 = 1.0047517553493037 + 6.300387486754831*mjd_current

theta = wrap_to_2pi(θ2)


# epc = Epoch(2020,11,9,20,31,4,0,tsys = "UTC") - 86164.09890369035*20

# jd_current = jd(epc)

# θ = earth_rotation(epc)
# era = iauEra00(MJD_ZERO, mjd(epc, tsys="UTC"))

# @show theta-era


# function mul_w(wp1,wp2,x)
#     return wp1*x + wp2*x
# end

# era_0 = 0.7064625862798124
# mjd_0 = 59154
era_0 = 1.7557955403696752
mjd_0 = 59215


delta_mjd = mjd_current - mjd_0

w_earth_rad_per_day_p1 = 6.30038
w_earth_rad_per_day_p2 = .7486754831e-5

θ_3 = era_0 + 6.300387486754831*delta_mjd
# θ_3 = era_0 + w_earth_rad_per_day_p1*delta_mjd + w_earth_rad_per_day_p2*delta_mjd

theta2 = wrap_to_2pi(θ2)


# @show era
#
# era = iauEra00(MJD_ZERO, mjd(epc + 86400, tsys="UT1"))
# @show era
#
# era = iauEra00(MJD_ZERO, mjd(epc, tsys="UT1")+1)
# @show era
#
#
# era_0 =

# mjd1 = 59161
#
# mjd2 = 1.85491


# earth_period = 2*pi/OMEGA_EARTH


# x = 0:1:60000
# y = wrap_to_2pi.(6.300387486754831*x)
#
#
# potent = 52962
# mat"
# figure
# hold on
# plot($x,$y)
# plot($x,2*pi*ones(length($x),1))
# hold off
# "



# get earth rotation angle

# MJD_int = 59161
# MJD_float = 1.85491
#
#     era_0 = 0.7064625862798124
#     mjd_0 = 59154
#
#     delta_mjd_int = MJD_int - mjd_0
#     delta_mjd_float = MJD_float
#
#     w_earth_rad_per_day_p1 = 6.30038
#     w_earth_rad_per_day_p2 = .7486754831e-5
#
#     a = delta_mjd_int*w_earth_rad_per_day_p1 + delta_mjd_int*w_earth_rad_per_day_p2
#     b = delta_mjd_float*w_earth_rad_per_day_p1 + delta_mjd_float*w_earth_rad_per_day_p2 + era_0
#
#     # here I get sin and cosine of GMST
#     sin_theta = sin(a)*cos(b) + cos(a)*sin(b)
#     cos_theta = cos(a)*cos(b) - sin(a)*sin(b)
#
# @show sin_theta
# @show cos_theta


function rotz(θ)
    """Rotation matrix for rotation about the z axis"""
    return [ cos(θ)  sin(θ) 0
            -sin(θ)  cos(θ) 0
             0       0      1]
end
function myECItoECEF(Era)
    """ECI to ECEF rotation matrix.

    Args:
        seconds_from_9_1_2020_to_epoch: (s)
        seconds_from_epoch: (s) will be time.monotonic()

    Returns:
        ECEF_Q_ECI: DCM
    """
    # t = seconds_from_9_1_2020_to_epoch + seconds_from_epoch
    return rotz(Era)
end

function rvecef_from_eci(r_eci,v_eci,Era)

    ecef_Q_eci = myECItoECEF(Era)

    r_ecef = ecef_Q_eci * r_eci

    ω = [0;0;OMEGA_EARTH]

    v_ecef = ecef_Q_eci*v_eci - ω × r_ecef

    return r_ecef, v_ecef
end

function rveci_from_ecef(r_ecef,v_ecef,Era)

    ecef_Q_eci = myECItoECEF(Era)
    eci_Q_ecef = transpose(ecef_Q_eci)

    r_eci = eci_Q_ecef*r_ecef

    ω = [0;0;OMEGA_EARTH]

    v_eci = eci_Q_ecef*v_ecef + ω × r_eci

    return r_eci, v_eci
end


oe0  = [R_EARTH + 400*1e3, 0.01,
            23, 0, 0, 34]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

r_eci = eci0[1:3]
v_eci = eci0[4:6]

r_ecef1,v_ecef1 = rvecef_from_eci(r_eci,v_eci,theta)

r_eci1,v_eci1 = rveci_from_ecef(r_ecef1,v_ecef1,theta)

@test isapprox(r_eci,r_eci1,rtol = 1e-6)
@test isapprox(v_eci, v_eci1, rtol = 1e-6)
