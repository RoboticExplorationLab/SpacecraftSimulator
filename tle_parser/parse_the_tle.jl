# using SatelliteToolbox

cd("/Users/kevintracy/devel/SpacecraftSimulator/tle_parser/sat_toolkit_env")
Pkg.activate(".")
using SatelliteToolbox
# read the TLE's
# line1 = "1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991"
# line2 = "2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102"

#ISS
line1 = "1 25544U 98067A   20259.77320602  .00000289  00000-0  13397-4 0  9991"
line2 = "2 25544  51.6439 260.1962 0000953 103.1904 238.6746 15.48949468246052"

TLEs = read_tle_from_string(line1, line2)

# initialize propagator
orbp = init_orbit_propagator(Val(:sgp4), TLEs[1])

sgp4d = orbp.sgp4d

(r,v) = sgp4!(sgp4d,0.0)
# o,r,v = propagate!(orbp,0.0)
