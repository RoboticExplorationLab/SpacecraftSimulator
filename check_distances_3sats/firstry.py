from sgp4.api import Satrec
import numpy as np

# 47463
one_l1 = "1 47463U 21006BC  21035.92088741  .00000935  00000-0  59038-4 0  9999"
one_l2 = "2 47463  97.5022  99.3873 0008943 212.6028 147.4651 15.11435868  2018"

#47467
two_l1 = "1 47467U 21006BG  21035.92039595  .00000941  00000-0  59318-4 0  9994"
two_l2 = "2 47467  97.5027  99.3886 0009009 210.7718 149.2985 15.11501890  2019"

# 47524
three_l1 = "1 47524U 21006DQ  21035.91983439  .00000936  00000-0  58917-4 0  9990"
three_l2 = "2 47524  97.5038  99.3910 0008991 208.4725 151.6015 15.11577583  2017"

# let me test with my own stuff
one_satellite = Satrec.twoline2rv(one_l1, one_l2)
two_satellite = Satrec.twoline2rv(two_l1, two_l2)
three_satellite = Satrec.twoline2rv(three_l1, three_l2)

# feb 1 2021 00:00:00
jd_p1 = 2459246.50000
jd_p2 = 0.0

# use sgp4 to get the current r_eci and v_eci (units of km and km/s)
e, r1, v1 = one_satellite.sgp4(jd_p1, jd_p2)
e, r2, v2 = two_satellite.sgp4(jd_p1, jd_p2)
e, r3, v3 = three_satellite.sgp4(jd_p1, jd_p2)

r1 = np.asarray(r1)
r2 = np.asarray(r2)
r3 = np.asarray(r3)
print(np.linalg.norm(r1 - r2))
print(np.linalg.norm(r2 - r3))
print(np.linalg.norm(r3 - r1))
