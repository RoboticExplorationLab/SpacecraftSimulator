# import brahe
from sgp4.api import Satrec
import time
import pysofa2
import pdb
import math
import numpy as np


def jd_from_tle(line1):
    """Take the year and day from the TLE and convert it to the two part JD.

    Args:
        line1: line 1 of the TLE (string)

    Returns:
        jd_p1: first part of jd
        jd_p2: second part of jd

    Comments: jd = jd_p1 + jd_p2
    """

    # get epoch info from TLE
    year = float(line1[18:20])
    day = float(line1[20:32])

    # find jd for given year
    if int(year) == 20:
        jd_p1 = 2.4588495e6

    if int(year) == 21:
        jd_p1 = 2.4592155e6

    # add the day (adjust for year starting on day 1 not day 0)
    jd_p2 = day - 1

    return jd_p1, jd_p2


def rvecef_from_eci(r_eci, v_eci, era):
    """Get rv in ecef from eci and earth rotation angle.
    Args:
        r_eci: position in eci (km)
        v_eci: velocity in eci wrt eci (km/s)
        era: earth rotation angle (radians)
    Returns:
        r_ecef: position in ecef (km)
        v_ecef: velocity in ecef wrt ecef (km/s)
    """

    sin_theta = math.sin(era)
    cos_theta = math.cos(era)

    r_ecef = np.array([
        cos_theta * r_eci[0] + sin_theta * r_eci[1],
        -sin_theta * r_eci[0] + cos_theta * r_eci[1],
        r_eci[2],
    ])

    omega_earth = 7.292115146706979e-5
    v_ecef = np.array([
        cos_theta * v_eci[0] + sin_theta * v_eci[1] + omega_earth * r_ecef[1],
        -sin_theta * v_eci[0] + cos_theta * v_eci[1] - omega_earth * r_ecef[0],
        v_eci[2],
    ])
    v_ecef = np.array([
        cos_theta * v_eci[0] + sin_theta * v_eci[1],
        -sin_theta * v_eci[0] + cos_theta * v_eci[1],
        v_eci[2],
    ]) + np.cross(-np.array([0, 0, omega_earth]), r_ecef)

    return r_ecef, v_ecef


def rveci_from_ecef(r_ecef, v_ecef, era):
    """Get rv in eci from ecef and earth rotation angle.
    Args:
        r_ecef: position in ecef (km)
        v_ecef: velocity in ecef wrt ecef (km/s)
        era: earth rotation angle (radians)
    Returns:
        r_eci: position in eci (km)
        v_eci: velocity in eci wrt eci (km/s)
    """

    sin_theta = math.sin(era)
    cos_theta = math.cos(era)

    r_eci = np.array([
        cos_theta * r_ecef[0] - sin_theta * r_ecef[1],
        sin_theta * r_ecef[0] + cos_theta * r_ecef[1],
        r_ecef[2],
    ])

    omega_earth = 7.292115146706979e-5
    v_eci = np.array([
        cos_theta * v_ecef[0] - sin_theta * v_ecef[1] - omega_earth * r_eci[1],
        sin_theta * v_ecef[0] + cos_theta * v_ecef[1] + omega_earth * r_eci[0],
        v_ecef[2],
    ])

    return r_eci, v_eci


# TLE's
line1 = "1 44400U 19038Q   20311.60590416  .00001180  00000-0  71834-4 0  9998"
line2 = "2 44400  97.5565 271.1792 0025309   0.8441 359.2825 15.12118200 74057"
satellite = Satrec.twoline2rv(line1, line2)

# get JD from TLE
jd_p1, jd_p2 = jd_from_tle(line1)

# full jd at time of TLE epoch
jd_tle = jd_p1 + jd_p2

print("jd at epoch", jd_tle)
# start sgp4 to get r and v
e, r_eci, v_eci = satellite.sgp4(jd_p1, jd_p2)

print("r_eci:", r_eci)
print("v_eci:", v_eci)

# get earth rotation angle from C SOFA library
ERA = pysofa2.Era00(jd_p1, jd_p2)
print("earth rotation angle", ERA)

# ecef
r_ecef, v_ecef = rvecef_from_eci(r_eci, v_eci, ERA)

print("r_ecef:", r_ecef)
print("v_ecef:", v_ecef)

# get current time
current_time = time.gmtime()

# pull JD from this time using C SOFA library
jd_current_p1, jd_current_p2 = pysofa2.Dtf2d(
    "UTC",
    current_time.tm_year,
    current_time.tm_mon,
    current_time.tm_mday,
    current_time.tm_hour,
    current_time.tm_min,
    current_time.tm_sec,
)

# calculate time since TLE epoch using julian day = 86400 seconds
seconds_since_epoch = 86400 * (jd_current_p1 + jd_current_p2 - jd_tle)

print("Seconds since TLE epoch", seconds_since_epoch)

# pdb.set_trace()

# # -- testing --
# R_EARTH = 6371
# ERA = 0  # random
# r_ecef = np.array([1e6, 0, 0])
# v_ecef = np.array([0, 0, 0])
# r_eci, v_eci = rveci_from_ecef(r_ecef, v_ecef, ERA)
# r_ecef2, v_ecef2 = rvecef_from_eci(r_eci, v_eci, ERA)
# print(r_ecef - r_ecef2)
# print(v_ecef - v_ecef2)

# R_EARTH = 6371
# ERA = 3.45  # random
# r_ecef = np.array([1e6, 2e6, 7.4e6])
# v_ecef = np.array([4, 3, 2.1])
# r_eci, v_eci = rveci_from_ecef(r_ecef, v_ecef, ERA)
# r_ecef2, v_ecef2 = rvecef_from_eci(r_eci, v_eci, ERA)
# print(r_ecef - r_ecef2)
# print(v_ecef - v_ecef2)

# R_EARTH = 6371
# ERA = 4.87  # random
# r_ecef = np.array([5e6, 7e6, -3e7])
# v_ecef = np.array([3, 4, -5.9])
# r_eci, v_eci = rveci_from_ecef(r_ecef, v_ecef, ERA)
# r_ecef2, v_ecef2 = rvecef_from_eci(r_eci, v_eci, ERA)
# print(r_ecef - r_ecef2)
# print(v_ecef - v_ecef2)
