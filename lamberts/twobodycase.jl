using SatelliteDynamics, ForwardDiff, MATLAB, LinearAlgebra
using StaticArrays


# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 5000e3, 0.3, 75.0, 45.0, 30.0, 0.0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

# Set the propagation end time to one orbit period after the start
T    = orbit_period(oe0[1])
epcf = epc0 + T

# Initialize State Vector
orb  = EarthInertialState(epc0, eci0, dt=1.0,
            mass=1.0, n_grav=0, m_grav=0,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
)

# Propagate the orbit
t, epc, eci_hist = sim!(orb, epcf)


r1 = eci_hist[1:3,1]
t1 = 0
r2 = eci_hist[1:3,2000]
t2 = t[2000]

# mat"
# [x,y,z] = sphere(20);
# figure
# hold on
# plot3($eci_hist(1,:),$eci_hist(2,:),$eci_hist(3,:),'linewidth',5)
# plot3($r1(1),$r1(2),$r1(3),'r*')
# plot3($r2(1),$r2(2),$r2(3),'r*')
# surf($R_EARTH*x,$R_EARTH*y,$R_EARTH*z)
# hold off
# "


function rk4_orbital(f::Function, t_n, x_n, u, h::Real)
    """Runge-Kutta 4th order integration. Epoch for time.
    Args:
        ODE:           f(t,x,u)        | function
        initial time:  t_n             | epoch or float64
        initial cond:  x_n             | vec
        control input: u               | vec
        time step:     h               | scalar
    Returns:
        x_{n+1}
    """

    k1 = h*f(t_n,x_n,u)
    k2 = h*f(t_n+h/2,x_n+k1/2,u)
    k3 = h*f(t_n+h/2,x_n+k2/2,u)
    k4 = h*f(t_n+h,x_n+k3,u)

    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end

function dynamics(t,x,u)
    r = SA[x[1],x[2],x[3]]
    a = -GM_EARTH*r/(norm(r)^3)
    return SA[x[4],x[5],x[6],a[1],a[2],a[3]]
end


dt = 1
N = 2000
X = [@SVector zeros(6) for i = 1:N ]

function residual(v)
    # X[1] = SVector{6}([r1;v])
    x =  SVector{6}([r1;v])
    for i = 1:(N-1)
        x = rk4_orbital(dynamics,0,x,0,dt)
    end
    return (x[1:3]-r2)
end

v1 = eci0[4:6]

# @btime residual(v1)

# JJ =  FiniteDiff.finite_difference_jacobian(residual,v1+100*randn(3))

function GN(v0)

# v0 = v1 + 100*randn(3)

new_S = 1e15
for i = 1:5

    res = residual(v0)
    S = dot(res,res)
    # @show S

    J = FiniteDiff.finite_difference_jacobian(residual,v0)
    newton_step = -J\res
    α = 1.0
    for i = 1:10
        new_v = v0 + α*newton_step
        new_res = residual(new_v)
        new_S = dot(new_res,new_res)
        if new_S < S
            v0 = copy(new_v)
            break
        else
            α /=2
        end
    end
    res = residual(v0)
    S = dot(res,res)
    @show S
    if norm(newton_step)<0.1
        break
    end

end

end
# v0 = v1 + 100*randn(3)
h = normalize(cross(r1,r2))
vcirc = sqrt(GM_EARTH/norm(r1))
v0 = normalize(cross(h,r1))*vcirc

true_v = eci_hist[4:6,2000]
GN(v0)



@show "yes"
# function lambert_min_energy(r0,r)
#
# tm = 1
#
#
# rn = norm(r)
# r0n = norm(r0)
# cosΔν = dot(r0,r)/(rn*r0n)
# sinΔν = tm*sqrt(1 - cosΔν )
# Δν = atan(sinΔν,cosΔν)
#
# c = sqrt(r0n^2 + rn^2 - 2*r0n*rn*cosΔν)
#
# s = (rn + r0n + c)/2
#
# a_min = s/2
# p_min = ((rn*r0n)/c)*(1 - cosΔν )
#
# e_min = sqrt(1 - (2*p_min)/(s))
#
# α_e = pi
#
# sinbe2 = sqrt((s-c)/s)
#
# # tmin_amin = sqrt(a_min^3/GM_EARTH)*()
#
# @show a_min
# @show e_min
#
# end
# lambert_min_energy(r1,r2)
