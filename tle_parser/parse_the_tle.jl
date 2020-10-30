# using SatelliteToolbox

cd("/Users/kevintracy/devel/SpacecraftSimulator/tle_parser/sat_toolkit_env")
Pkg.activate(".")
using SatelliteToolbox, SatelliteDynamics, HTTP
const SD = SatelliteDynamics
# read the TLE's
# line1 = "1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991"
# line2 = "2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102"

#ISS
# line1 = "1 25544U 98067A   20267.41567050  .00001718  00000-0  39590-4 0  9992"
# line2 = "2 25544  51.6436 222.3966 0001436  98.6522  47.1481 15.48789843247244"

r = HTTP.get("https://www.celestrak.com/NORAD/elements/cubesat.txt")
ss = String(r.body)

name_range = findfirst("SONATE",ss)

id = name_range[1]

line1 = ss[id+26:id+94]
line2 = ss[id+97:id+97+68]
# SONATE NORAD = 44400
# line1 = "1 44400U 19038Q   20267.59883466  .00000414  00000-0  27750-4 0  9993"
# line2 = "2 44400  97.5499 227.6445 0025856 145.1208 215.1722 15.12040781 67400"
# line1 =

TLEs = read_tle_from_string(line1, line2)

# initialize propagator
orbp = init_orbit_propagator(Val(:sgp4), TLEs[1])

sgp4d = orbp.sgp4d

# this gives us the R and V in eci
(r_tle,v_tle) = sgp4!(sgp4d,0.0)
# o,r,v = propagate!(orbp,0.0)
jd_epoch = sgp4d.epoch
# mjd_epoch = jd_epoch - MJD_ZERO
# mjd_epoch = jd_epoch - 2.4000005e6

date_tuple = SD.jd_to_caldate(jd_epoch )

@show r_tle*1000
@show v_tle*1000
@show jd_epoch

# T = Epoch()
