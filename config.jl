using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox

const SD = SatelliteDynamics
const ST = SatelliteToolbox

J = [1.959e-4 2016.333e-9 269.176e-9;
    2016.333e-9 1.999e-4 2318.659e-9;
    269.176e-9 2318.659e-9 1.064e-4]

sc_mass = 177.808e-3



line1 = "1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991"
line2 = "2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102"


ST.read_tle_from_string(line1, line2)
