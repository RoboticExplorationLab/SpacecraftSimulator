import numpy as np

import time


class igrfclass:
    gh = [
        -29404.8,
        -1450.9,
        4652.5,
        -2499.6,
        2982.0,
        -2991.6,  # 2020
        1677.0,
        -734.6,
        1363.2,
        -2381.2,
        -82.1,
        1236.2,  # 2020
        241.9,
        525.7,
        -543.4,
        903.0,
        809.5,
        281.9,  # 2020
        86.3,
        -158.4,
        -309.4,
        199.7,
        48.0,
        -349.7,  # 2020
        -234.3,
        363.2,
        47.7,
        187.8,
        208.3,
        -140.7,  # 2020
        -121.2,
        -151.2,
        32.3,
        13.5,
        98.9,
        66.0,  # 2020
        65.5,
        -19.1,
        72.9,
        25.1,
        -121.5,
        52.8,  # 2020
        -36.2,
        -64.5,
        13.5,
        8.9,
        -64.7,
        68.1,  # 2020
        80.6,
        -76.7,
        -51.5,
        -8.2,
        -16.9,
        56.5,  # 2020
        2.2,
        15.8,
        23.5,
        6.4,
        -2.2,
        -7.2,  # 2020
        -27.2,
        9.8,
        -1.8,
        23.7,
        9.7,
        8.4,  # 2020
        -17.6,
        -15.3,
        -0.5,
        12.8,
    ]

    def __init__(self):
        self.NN = 33

    def testfx(self, N):
        print(self.gh)
        print(N)


igg = igrfclass()

print(igg.gh)

igg.testfx(69)


# a = np.zeros(1000)
# t1 = time.time()
# for i in range(10000):
#     a *= 0

# print(time.time() - t1)

# a = np.zeros(1000)
# t1 = time.time()
# for i in range(10000):
#     a = a * 0

# print(time.time() - t1)
class Earth:
    mu = 3.986004418e5  # specific gravitational parameter
    R = 6371.009  # radius in km
    J2 = 1.08262668e-3  # J2 value
    e = 0.08181919  # elliptical Earth


earth = Earth()

print(earth.mu)


import math 

ground_stations = [[37.4241,-122.166],[37.4241,-75],[-30.5595,22.9375],[90,0]]

class Earth():
	"""
	class for storing all Earth-related parameters
	"""
	mu = 3.986004418E5 # specific gravitational parameter
	R = 6371.009 # radius in km
	J2 = 1.08262668E-3 # J2 value
	e = 0.08181919 # elliptical Earth


earth = Earth()

for station in ground_stations:
	# constants
	deg2rad = math.pi/180
	# get station ECEF
	N = earth.R / math.sqrt(1 - earth.e ** 2 * math.sin(station[0]*deg2rad) ** 2)
	r_station = np.array(
	[N * np.vector.cos(station[0]*deg2rad) * np.vector.cos(station[1]*deg2rad),
	N * np.vector.cos(station[0]*deg2rad) * np.vector.sin(station[1]*deg2rad),
	N*(1-earth.e ** 2) * np.vector.sin(station[0]*deg2rad)])
    print(r_station)