import math
import pdb
import numpy as np


def mjd_from_gps(GNSS_week, TOW):
    """Returns MJD from gps time.

    Args:
        GNSS_week: Weeks since 0h January 6th, 1980 (uint16, units: weeks)
        TOW: Seconds into the week (uint16, units: 1/100 seconds)

    Returns:
        MJD_current: MJD as described by the GPS time

    Comments:
        Both of the inputs to this function are the raw GPS time parameters.
        The TOW term is scaled by 0.01 after being input to the function.

        This can be tested with <http://leapsecond.com/java/gpsclock.htm>
    """

    # seconds since the week started (GPS is off from UTC by 18 seconds)
    TOW = (TOW // 100) - 18

    # MJD when GPS time started
    MJD_gps_epoch = 44244

    # get days since this epoch
    GNSS_days = GNSS_week * 7 + (TOW / 86400)

    # get current MJD
    MJD_current = MJD_gps_epoch + GNSS_days

    return MJD_current


def gps_from_mjd(MJD_current):
    """Returns GPS time from MJD.

    Args:
        MJD_current: MJD as described by the GPS time
        GNSS_week: Weeks since 0h January 6th, 1980 (uint16, units: weeks)
        TOW: Seconds into the week (uint16, units: 1/100 seconds)

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


def test_gps_time():
    """Five tests run from http://leapsecond.com/java/gpsclock.htm"""

    # -------------- test 1 -----------------
    # weeks from midnight january 6th 1980
    GNSS_week = 2131

    # seconds since the week started
    TOW = 15789200

    # given MJD
    MJD = 59162.8274

    mjd_current = mjd_from_gps(GNSS_week, TOW)
    np.testing.assert_approx_equal(mjd_current, MJD, 8)
    GNSS_week2, TOW2 = gps_from_mjd(mjd_current)

    np.testing.assert_equal(GNSS_week, GNSS_week2)
    np.testing.assert_equal(TOW, TOW2)

    # -------------- test 2 -----------------

    # weeks from midnight january 6th 1980
    GNSS_week = 2131

    # seconds since the week started
    TOW = 16028200

    # given MJD
    MJD = 59162.8549

    mjd_current = mjd_from_gps(GNSS_week, TOW)
    np.testing.assert_approx_equal(mjd_current, MJD, 8)
    GNSS_week2, TOW2 = gps_from_mjd(mjd_current)

    np.testing.assert_equal(GNSS_week, GNSS_week2)
    np.testing.assert_equal(TOW, TOW2)

    # -------------- test 3 -----------------

    # weeks from midnight january 6th 1980
    GNSS_week = 2131

    # seconds since the week started
    TOW = 19047200

    # given MJD
    MJD = 59163.20432

    mjd_current = mjd_from_gps(GNSS_week, TOW)
    np.testing.assert_approx_equal(mjd_current, MJD, 8)
    GNSS_week2, TOW2 = gps_from_mjd(mjd_current)

    np.testing.assert_equal(GNSS_week, GNSS_week2)
    np.testing.assert_equal(TOW, TOW2)

    # -------------- test 4 -----------------

    # weeks from midnight january 6th 1980
    GNSS_week = 2131

    # seconds since the week started
    TOW = 41970300

    # given MJD
    MJD = 59165.85746

    mjd_current = mjd_from_gps(GNSS_week, TOW)
    np.testing.assert_approx_equal(mjd_current, MJD, 8)
    GNSS_week2, TOW2 = gps_from_mjd(mjd_current)

    np.testing.assert_equal(GNSS_week, GNSS_week2)
    np.testing.assert_equal(TOW, TOW2)

    # -------------- test 5 -----------------

    # weeks from midnight january 6th 1980
    GNSS_week = 2132

    # seconds since the week started
    TOW = 51638600

    # given MJD
    MJD = 59173.97648

    mjd_current = mjd_from_gps(GNSS_week, TOW)
    np.testing.assert_approx_equal(mjd_current, MJD, 8)
    GNSS_week2, TOW2 = gps_from_mjd(mjd_current)

    np.testing.assert_equal(GNSS_week, GNSS_week2)
    np.testing.assert_equal(TOW, TOW2)


test_gps_time()

print("yes")
