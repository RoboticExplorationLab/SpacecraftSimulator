from sgp4.api import Satrec
import time
import pysofa2
import pdb
import math
import numpy as np
import string


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

    return r_ecef, v_ecef


def gps_from_mjd(MJD_current):
    """Returns GPS time from MJD.

    Args:
        MJD_current: MJD as described by the GPS time

    Returns:
        GNSS_week: Weeks since 0h January 6th, 1980 (uint16, units: weeks)
        TOW: Seconds into the week (uint16, units: 1/100 seconds)

    Comments:
        Both of the outputs to this function are the raw GPS time parameters.
        The TOW term is scaled by 0.01 on the way out of the function.

        This can be tested with <http://leapsecond.com/java/gpsclock.htm>
    """

    # MJD when GPS time started
    MJD_gps_epoch = 44244

    # julian days since this epoch
    GNSS_days = MJD_current - MJD_gps_epoch

    # number of weeks since this epoch
    GNSS_week_float = GNSS_days / 7

    # convert this week to the floor int
    GNSS_week = int(np.floor(GNSS_week_float))

    # convert the excess to seconds, and add 18 seconds for current UTC offset
    TOW = (GNSS_week_float % 1) * 7 * 86400 + 18

    # convert to the scaling (0.01s) that GPS uses
    TOW = round(TOW * 100)

    return GNSS_week, TOW


# TLE's (Sonate from <https://www.celestrak.com/NORAD/elements/cubesat.txt>)
line1 = "1 44400U 19038Q   20345.94933858  .00001475  00000-0  88741-4 0  9995"
line2 = "2 44400  97.5639 305.1881 0023394 237.5108 122.3864 15.12211443 79248"

# line1 = "1 99789U          21022.64108674  .00000009  00000-0  59233-6 0 00008"
# line2 = "2 99789 097.4984 077.3576 0009198 294.5224 091.6031 15.08170883000346"


def checksum(line):
    # L = string.strip(line)
    L = line1
    cksum = 0
    for i in range(68):
        c = L[i]
        if c == ' ' or c == '.' or c == '+' or c in string.ascii_letters:
            continue
        elif c == '-':
            cksum = cksum + 1
        else:
            cksum = cksum + int(c)

    cksum %= 10

    return cksum


checksum(line1)
pdb.set_trace()

# let me test with my own stuff
satellite = Satrec.twoline2rv(line1, line2)

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

jd_current_p1 = 2.459237141086736e6
jd_current_p2 = 0.0
# use sgp4 to get the current r_eci and v_eci (units of km and km/s)
sgp4_t1 = time.time()
e, r_eci, v_eci = satellite.sgp4(jd_current_p1, jd_current_p2)
print("time for propagation:", time.time() - sgp4_t1)

print("---------------Propagator Data---------------")
print("r_eci (km):", r_eci)
print("v_eci (km/s):", v_eci)

# get earth rotation angle from C SOFA library
ERA = pysofa2.Era00(jd_current_p1, jd_current_p2)
print("earth rotation angle (radians)", ERA)

print("------------------GPS DATA-------------------")
# get ecef position and velocity
r_ecef, v_ecef = rvecef_from_eci(r_eci, v_eci, ERA)

print("r_ecef (km):", r_ecef)
print("v_ecef (km/s):", v_ecef)

# get GPS readings (same as raw sensor)
MJD_ZERO = 2.4000005e6
MJD_current = (jd_current_p1 + jd_current_p2) - MJD_ZERO
GNSS_week, TOW = gps_from_mjd(MJD_current)

print("GNSS Week (weeks):", GNSS_week)
print("TOW (0.01 seconds):", TOW)
