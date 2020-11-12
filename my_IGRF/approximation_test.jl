# function IGRFecef(r_ecef,epc)
#
#     # # get decimal date
#     date = SD.caldate(epc)
#     year = date[1]
#     day = SD.day_of_year(epc)
#     decimal_date = year + day/365.2425
#     #
#     # # eci and ecef location
#     # ecef_Q_eci = SD.rECItoECEF(epc)
#     # eci_Q_ecef = transpose(ecef_Q_eci)
#     #
#     # # get ecef location
#     # r_ecef = ecef_Q_eci*r_eci
#
#     # long lat geod
#     longitude,latitude,altitude = SD.sECEFtoGEOC(r_ecef,use_degrees=false)
#
#     # IGRF
#     # SatelliteToolbox v0.7.1
#     # B_ned_nT = igrf(decimal_date, norm(r_ecef), latitude, longitude, Val(:geocentric))
#     # SatelliteToolbox v0.6.3
#     B_ned_nT = igrf12(decimal_date, norm(r_ecef), latitude, longitude, Val{:geocentric})
#
#     # NED and ECEF stuff
#     ecef_Q_ned = ecef_Q_ned_mat(longitude,latitude)
#
#     # conver to eci
#     B_ecef_nT = ecef_Q_ned*B_ned_nT
#
#     # convert from nT to T
#     return B_ecef_nT*1e-9
# end

function spline(rv0,rv1,t)

    D = rv0[1:3]
    C = rv0[4:6]


    coeffs = (1/t^3)*[-2 t; 3*t -t^2]*[(rv1[1:3] - C*t - D)';
                                       (rv1[4:6] - C)'       ]

    A = vec(coeffs[1,:])
    B = vec(coeffs[2,:])

    return A,B,C,D
end


# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550e3, 0.0, 48, 20, 20, 0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

function fode(x)

    r = x[1:3]

    return [x[4:6];-GM_EARTH*r/norm(r)^3]
end


function rk4(f,x_n,dt)
k1 = dt*f(x_n)
k2 = dt*f(x_n+k1/2)
k3 = dt*f(x_n+k2/2)
k4 = dt*f(x_n+k3)

return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end


dt = 10
N = 30
X = fill(zeros(6),N)
X[1]=eci0

for i = 1:N-1

    X[i+1] = rk4(fode,X[i],dt)

end

ecef0 = SD.sECItoECEF(epc0,X[1])
ecef1 = SD.sECItoECEF(epc0+N*dt,X[end])

# linap = (-ecef0 + ecef1)/N

# MATB = []
#true igrf history
IGRF_true = fill(zeros(3),N)
IGRF_approx = fill(zeros(3),N)
# for i = 1:N
#
#     IGRF_true = IGRF13(X[i][1:3],epc+(i-1)*dt)
#     IGRF_fake = IGRFecef(ecef0 + (1-N)*linap,epc)

t = (N-1)*dt

A,B,C,D = spline(X[1],X[end],t)

function splinefx(A,B,C,D,t)
    return A*t^3 + B*t^2 + C*t + D
end

Xspl = copy(X)
for i = 1:N
    Xspl[i] = splinefx(A,B,C,D,(i-1)*dt)
end

X = mat_from_vec(X)
Xspl = mat_from_vec(Xspl)

er = X[1:3,:] - Xspl

e = zeros(size(er,2))

for i = 1:size(er,2)
    e[i] = norm(er[:,i])
end


function rad_from_arcsec(arcsec)
    return deg2rad(arcsec/3600)
end


# function R_x(theta)
#     return active_rotation_axis_angle([1.0;0;0],theta)
# end
# function R_y(theta)
#     return active_rotation_axis_angle([0;1.0;0],theta)
# end
# function R_z(theta)
#     return active_rotation_axis_angle([0;0;1.0],theta)
# end

# function ECEF_Q_ECI_mat(epc)
#
#     # convert to julian date
#     JD = jd(epc)
#
#     # julian centuries
#     T = (JD - 2451545.0)/36525.0
#
#     # three angles
#     ζ = 2306.2181*T  + 0.30188*T^2 + 0.017998*T^3
#
#     ν = 2004.3109*T - 0.42665*T^2 - 0.041833*T^3
#
#     z = ζ + .79280 *T^2 + 0.000205*T^3
#
#     # DCM
#     ECEF_Q_ECI = SD.Rz(-rad_from_arcsec(z))*SD.Ry(rad_from_arcsec(ν))*SD.Rz(-rad_from_arcsec(ζ))
#
#     return ECEF_Q_ECI
#
# end
#
jdut1 = jd(epc0)
#
tut1= ( jdut1 - 2451545.0 ) / 36525.0;
#
temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841
#
@show SD.Rz(deg2rad(temp))

# epc0 = Epoch(2000, 1, 1, 0, 0, 0, 0.0)


d = jd(epc0)-2451545.0

gmst = 280.4606 + 360.985643*d

SD.Rz(deg2rad(gmst))

Q_real = SD.rECItoECEF(epc0)
Q = SD.Rz(deg2rad(gmst))



N = 10000
er = zeros(N)

for i = 1:N
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

    epc0 += 100000*randn()
    d = jd(epc0)-2451545.0

    gmst = 280.4606 + 360.985643*d

    # SD.Rz(deg2rad(gmst))

    Q_real = SD.rECItoECEF(epc0)
    Q = SD.Rz(deg2rad(gmst))

r = normalize(randn(3))*(R_EARTH + 550e3)

er[i] = norm(Q*r - Q_real*r)

end



function K_size(nx,nu,N)
    return nx*nu*(N-1)
end

function x_size(nx,N)
    return 2*nx*N
end
function u_size(nu,N)
    return 2*nu*N
end

function siz(nx,nu,N)
    return K_size(nx,nu,N) + x_size(nx,N) +u_size(nu,N)
end

nx = 2
nu = 1
N = 450

siz(nx,nu,N)

nx = 6
nu = 3
N = 25
siz(nx,nu,N)
