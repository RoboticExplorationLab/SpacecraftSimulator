ss_sim_path =  dirname(dirname(dirname(@__FILE__)))

cd(joinpath(ss_sim_path,"virtual_env"))
Pkg.activate(".")

using SOFA

function rotz(θ)
    """Rotation matrix for rotation about the z axis"""
    return [ cos(θ)  sin(θ) 0
            -sin(θ)  cos(θ) 0
             0       0      1]
end
function myECItoECEF(delta_earth_rotation,seconds_from_epoch)
    """ECI to ECEF rotation matrix.

    Args:
        seconds_from_9_1_2020_to_epoch: (s)
        seconds_from_epoch: (s) will be time.monotonic()

    Returns:
        ECEF_Q_ECI: DCM
    """
    # t = seconds_from_9_1_2020_to_epoch + seconds_from_epoch
    return rotz(delta_earth_rotation + seconds_from_epoch*7.292115146706979e-5)
end

function mat_from_vec(a::Union{Array{Array{Float64,1},1},Array{Array{Float32,1},1}})
    "Turn a vector of vectors into a matrix"


    rows = length(a[1])
    columns = length(a)
    A = zeros(rows,columns)

    for i = 1:columns
        A[:,i] = a[i]
    end

    return A
end

function get_el(r_ecef,r_stanford)
    # returns degrees
    return 90 - abs(rad2deg(acos(dot(normalize(r_ecef - r_stanford),normalize(r_stanford)))))
end

function visible(el,on_sight,epoch,min_el_deg)

    # if you don't already know, we are on sight
    if !on_sight && el > min_el_deg
        println("Pass open")
        @show epoch
        on_sight = true
        return on_sight

    # if you think we are on sight, and we aren't change it
    elseif on_sight && el < min_el_deg
        println("Pass closed")
        @show epoch
        on_sight = false
        return on_sight

    # if you aren't suprised by the result, change nothing
    else
        return on_sight
    end
end




using LinearAlgebra, SatelliteDynamics, MATLAB, ProgressMeter

using LinearAlgebra


function FODE(epc,eci_state,u)

    # unpack state
    r_eci = eci_state[1:3]
    v_eci = eci_state[4:6]


    # J2 only
    a_eci = FODE_J2(r_eci)


    return [v_eci;a_eci]
end

function FODE_J2(r_eci)

    # relavent parameters
    J2 = J2_EARTH
    # J2 = 0.0
    μ = GM_EARTH

    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  (-μ/r^3)*[x*(1 - Re_r_sqr*(five_z_sqr - 1));
                      y*(1 - Re_r_sqr*(five_z_sqr - 1));
                      z*(1 - Re_r_sqr*(five_z_sqr - 3))]

end


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

    x_np1 = (x_n + (1/6)*(k1+2*k2+2*k3 + k4))



    return x_np1
end

function runit(r,v,jd_epoch)

# x0 = vec([6928,0,0,0,-.871817082027,7.5348948722])*1000


# x0 = vec([-4.724512772111592e6, -5.064712686679697e6, 27.601182963844654,
# -714.9400927271462, 692.5574062503943, 7514.06327176573])
x0 = [r_tle;v_tle]*1000
# x0 = [-4689.932188;-5096.192909;10.925068;-.703127;0.689545;7.515902]*1000
dt = 10.0

stanford_ecef = sGEODtoECEF([-122.1697;37.4275;0.0],use_degrees = true)
# date = jd_to_caldate(2.4591159156705e6)
# jd_epoch = 2.45911543705526e6
mjd_epoch = jd_epoch - MJD_ZERO

delta_earth_rotation = iauEra00(MJD_ZERO,mjd_epoch)
@show delta_earth_rotation
date = jd_to_caldate(jd_epoch)
# 2.45911543705526e6
epoch = Epoch(date[1],date[2],date[3],date[4],date[5],date[6],date[7];tsys = "UTC")
# epoch_anchor = epoch
println(".")
println("Starting Time")
# @show epoch
tf = 40*6400.0

t_vec = 0:dt:tf

#N = length(t_vec)
N = 2001
X = fill(zeros(6),N)

X[1] = x0
min_el_deg = 0.0
pre_on_sight = false
on_sight = false
pass_timing = []
el_for_passes = []
passes = 0
for i = 1:N-1


    ECEF_Q_ECI = myECItoECEF(delta_earth_rotation,dt*(i-1))
    # check pass
    # r_ecef = rECItoECEF(epoch)*X[i][1:3]
    r_ecef = ECEF_Q_ECI*X[i][1:3]
    el = get_el(r_ecef,stanford_ecef)
    on_sight = visible(el,on_sight,epoch,min_el_deg)

    # keep track of passes
    if !pre_on_sight && on_sight
        push!(el_for_passes,[])
        push!(pass_timing,[])
        passes += 1
        push!(pass_timing[passes],epoch)
        # pass_iter = 0
    end
    if on_sight
        # push!(pass_timing[passes])
        push!(el_for_passes[passes],el)
    end
    if !on_sight && pre_on_sight
        push!(pass_timing[passes],epoch)
    end


    pre_on_sight = copy(on_sight)


    X[i+1] = rk4_orbital(FODE, t_vec[i], X[i], zeros(3), dt)



    epoch += dt
end


return X, el_for_passes, pass_timing

end

X, el_for_passes, pass_timing = runit(r_tle,v_tle,jd_epoch)

X_mat = mat_from_vec(X)

max_el_each_pass = zeros(length(el_for_passes))
for i = 1:length(el_for_passes)
    max_el_each_pass[i] = maximum(el_for_passes[i])
    push!(pass_timing[i],(pass_timing[i][2] - pass_timing[i][1])/60)
    push!(pass_timing[i],max_el_each_pass[i])
end


# mat"figure
# hold on
# plot3($X(1,:),$X(2,:),$X(3,:))
# hold off
# "
# ten_am_el = el_for_passes[4]
# mat"
# figure
# hold on
# plot($ten_am_el)
# hold off
# "


# stanford_ecef = sGEODtoECEF([-122.1697;37.4275;0.0],use_degrees = true)
