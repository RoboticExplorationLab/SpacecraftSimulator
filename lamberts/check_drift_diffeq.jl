using LinearAlgebra
using SatelliteDynamics
using MATLAB
using Random
using OrdinaryDiffEq
using ProgressMeter

# canonical units
const tscale = 806.0
const dscale = copy(R_EARTH)


r_ecef = [701408.178, 6878040.638, 106243.787]
v_ecef = [1477.693, -269.846, 7527.745]

# epc = Epoch("2021-01-24T15:59:24.718Z")
const epcc = Epoch("2021-01-24T15:59:24.718Z")

eci_Q_ecef = rECEFtoECI(epcc)

r_eci = eci_Q_ecef*r_ecef
v_eci = eci_Q_ecef * (v_ecef + [0;0;OMEGA_EARTH] × r_ecef )

days = 14
T = 24*3600*days










function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
             area_srp::Real=1.0, coef_srp::Real=1.8,
             n_grav::Integer=10, m_grav::Integer=10)
    """Accelerations for spacecraft in LEO, ForwardDiff friendly"""

    # Extract position and velocity
    r = x[1:3]
    v = x[4:6]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E  = earth_rotation(epc)
    W  = polar_motion(epc)
    R  = W * E * PN

    # Compute sun and moon position
    r_sun  = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration (eltype(x) makes this forward diff friendly)
    a = zeros(eltype(x), 3)

    # spherical harmonic gravity
    a += accel_gravity(x, R, n_grav, m_grav)

    # atmospheric drag
    ρ = density_harris_priester(epc,r)
    a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

    # SRP
    nu = eclipse_conical(x, r_sun)
    a += nu*accel_srp(x, r_sun, mass, area_srp, coef_srp)

    # third body sun
    a += accel_thirdbody_sun(x, r_sun)

    # third body moon
    a += accel_thirdbody_moon(x, r_moon)

    return a
end
function dynamics!(dx,x,p,t)
    """ODE for orbital motion"""
    r1 = x[1:3]*dscale
    v1 = x[4:6]*(dscale/tscale)
    r2 = x[7:9]*dscale
    v2 = x[10:12]*(dscale/tscale)
    r3 = x[13:15]*dscale
    v3 = x[16:18]*(dscale/tscale)

    # the current time is (epc + t*dscale)
    t_current = epcc+t*tscale

    dx .= [v1 / (dscale/tscale);
            accel_perturbations(t_current,[r1;v1]) / (dscale/tscale^2);
            v2 / (dscale/tscale);
            accel_perturbations(t_current,[r2;v2]) / (dscale/tscale^2);
            v3 / (dscale/tscale);
            accel_perturbations(t_current,[r3;v3]) / (dscale/tscale^2)]
end




mc = []
trials = 3
t1 = time()
@showprogress "mc run" for i = 1:500
    eci1 = [r_eci;v_eci]
    eci2 = eci1 + [zeros(3);(.1+.03*randn())*normalize(randn(3))]
    eci3 = eci1 + [zeros(3);(.1+.03*randn())*normalize(randn(3))]
    eci1[1:3] /= dscale
    eci1[4:6] /= (dscale/tscale)
    eci2[1:3] /= dscale
    eci2[4:6] /= (dscale/tscale)
    eci3[1:3] /= dscale
    eci3[4:6] /= (dscale/tscale)
    x0 = [eci1;eci2;eci3]
    tspan = (0.0,T/tscale)
    prob = ODEProblem(dynamics!,x0,tspan)
    sol = solve(prob,Vern7(),reltol = 1e-6)

    x_1 = [dscale*sol.u[i][1:3] for i = 1:length(sol.u)]
    x_2 = [dscale*sol.u[i][7:9] for i = 1:length(sol.u)]
    x_3 = [dscale*sol.u[i][13:15] for i = 1:length(sol.u)]

    r_dist1 = [norm(x_1[i] - x_2[i]) for i = 1:length(sol.t)]
    r_dist2 = [norm(x_1[i] - x_3[i]) for i = 1:length(sol.t)]
    r_dist3 = [norm(x_2[i] - x_3[i]) for i = 1:length(sol.t)]
    reld = [r_dist1,r_dist2,r_dist3]
    t_vec = copy(Array(sol.t))

    push!(mc,(reld = reld,t_vec = t_vec))
end
@show time() - t1

# mat"
# figure
# hold on
# plot($t_vec/(24*3600),$r_dist1/1000)
# plot($t_vec/(24*3600),$r_dist2/1000)
# plot($t_vec/(24*3600),$r_dist3/1000)
# xlabel('Time (days)')
# ylabel('Distance (km)')
# hold off
# "

using JLD2
@save "mc_drift.jld2" mc

max_distances = zeros(length(mc))
for i = 1:length(mc)
    max_distances[i]=maximum(maximum.(mc[i].reld))
end

mat"
figure
hold on
title('Maximum Drift Distance After 14 Days (500 Trials)')
histogram($max_distances/1000,20)
xlabel('Max Drift Distance (km)')
ylabel('Frequency')
saveas(gcf,'drift_mc.png')
"
