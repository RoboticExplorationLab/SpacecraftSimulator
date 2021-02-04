using LinearAlgebra
using SatelliteDynamics
using MATLAB

const tscale = 806.0
const dscale = copy(R_EARTH)
function run_sim(eci_initial,epc,epcf)
        orb  = EarthInertialState(epc, eci_initial, dt=30.0,
                mass=100.0, n_grav=20, m_grav=20,
                drag=true, srp=true,
                moon=true, sun=true,
                relativity=false
        )
        return sim!(orb, epcf)
end
r_ecef = [701408.178, 6878040.638, 106243.787]
v_ecef = [1477.693, -269.846, 7527.745]

# epc = Epoch("2021-01-24T15:59:24.718Z")
const epcc = Epoch("2021-01-24T15:59:24.718Z")

eci_Q_ecef = rECEFtoECI(epc)

r_eci = eci_Q_ecef*r_ecef

v_eci = eci_Q_ecef * (v_ecef + [0;0;OMEGA_EARTH] × r_ecef )

eci1 = [r_eci;v_eci]
eci2 = eci1 + [zeros(3);.1*normalize(randn(3))]
eci3 = eci1 + [zeros(3);.1*normalize(randn(3))]
days = 20
T = 24*3600*days
epcf = epc + T

# # Propagate the orbit
# t_hist1, epc_hist1, eci_hist1 = run_sim(eci1,epc,epcf)
# t_hist2, epc_hist2, eci_hist2 = run_sim(eci2,epc,epcf)
# t_hist3, epc_hist3, eci_hist3 = run_sim(eci3,epc,epcf)
#
#
r_dist1 = [norm(eci_hist1[1:3,i] - eci_hist2[1:3,i]) for i = 1:length(t_hist1)]
r_dist2 = [norm(eci_hist1[1:3,i] - eci_hist3[1:3,i]) for i = 1:length(t_hist1)]
r_dist3 = [norm(eci_hist2[1:3,i] - eci_hist3[1:3,i]) for i = 1:length(t_hist1)]
#
mat"
figure
hold on
plot($t_hist1/(24*3600),$r_dist1/1000)
plot($t_hist1/(24*3600),$r_dist2/1000)
plot($t_hist1/(24*3600),$r_dist3/1000)
xlabel('Time (days)')
ylabel('Distance (km)')
hold off
"


# function diffeq_ode(dx,x,p,t)
#         r1 = x[1:3]*dscale
#         v1 = x[4:6]*(dscale/tscale)
#
#         fderiv_earth_orbit(epcc + t, x::Array{<:Real} ;
#                      mass::Real=1.0, area_drag::Real=1.0, coef_drag::Real=2.3,
#                      area_srp::Real=1.0, coef_srp::Real=1.8,
#                      n_grav::Integer=20, m_grav::Integer=20,
#                      drag::Bool=true, srp::Bool=true, moon::Bool=true, sun::Bool=true,
#                      relativity::Bool=false)
