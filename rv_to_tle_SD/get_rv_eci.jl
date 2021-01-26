using SatelliteDynamics, LinearAlgebra, JLD2

epc = Epoch("2021-01-22T15:23:09.894Z")

r_ecef = [682334.037, 6880772.409, 12376.330]

v_ecef = [1486.455, -162.866, 7528.889]


eci_Q_ecef = rECEFtoECI(epc)

r_eci = eci_Q_ecef*r_ecef

v_eci = eci_Q_ecef * v_ecef + eci_Q_ecef*cross([0;0;OMEGA_EARTH],r_ecef )

jd_epc = jd(epc)

sc_info = (r_eci = r_eci,v_eci = v_eci, jd_epc = jd_epc)

eci0 = [r_eci;v_eci]

sCARTtoOSC(eci0,use_degrees = true)

# using JLD2
# @save "spacecraft_info.jld2" sc_info
#
#
# # TLE from STK
#
# 1 99999U          21022.64108674  .00000009  00000-0  59233-6 0 00008
# 2 99999 097.4984 077.3576 0009198 294.5224 091.6031 15.08170883000346
