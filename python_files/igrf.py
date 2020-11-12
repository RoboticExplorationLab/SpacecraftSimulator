import numpy as np
import math


rad2deg = 180 / math.pi
deg2rad = math.pi / 180
R_EARTH = 6.378136300e6  # m

date = 2020.3
latitude_degrees = 33
elongitude_degrees = 76
r_norm_km = 1.05 * R_EARTH / 1000


def igrf13_5(date, latitude_degrees, elongitude_degrees, r_norm_km):
    """IGRF 13 model truncated to order/degree 5.

    References:
        [1]: https://www.ngdc.noaa.gov/IAGA/vmod/igrf13.f

    Arguments:
        date: decimal date
        latitude_degrees: geocentric latitude in degrees
        elongitude_degrees: geocentric east longitude in degrees
        r_norm_km: radius from center of earth to spacecraft in km

    Outputs:
        b_ecef_T: magnetic field expressed in ECEF frame, in Teslas (T)
    """

    # IGRF13 Coefficients
    gh = np.array(
        [
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
    )

    # order = 5

    # colat = colatd_from_latd(lat)
    colat = 90 - latitude_degrees

    # Declaration of variables
    fn = 0
    gn = 0
    kmx = 0
    ll = 0
    nc = 0
    nmx = 0
    x = 0.0
    y = 0.0
    z = 0.0
    t = 0.0
    tc = 0.0

    t = date - 2020
    tc = 1.0

    ll = 0
    nmx = 5

    # nc = int(nmx * (nmx + 2))
    nc = 35

    # kmx = int((nmx + 1) * (nmx + 2) / 2)
    kmx = 21

    # allocate
    cl = np.zeros(nmx)
    sl = np.zeros(nmx)
    p = np.zeros(kmx)
    q = np.zeros(kmx)

    r = r_norm_km
    ct = math.cos(colat * math.pi / 180)
    st = math.sin(colat * math.pi / 180)
    cl[0] = math.cos(elongitude_degrees * math.pi / 180)
    sl[0] = math.sin(elongitude_degrees * math.pi / 180)
    Cd = 1.0
    sd = 0.0
    l = 1
    m = 1
    n = 0

    ratio = 6371.2 / r
    rr = ratio ** 2

    p[0] = 1.0
    p[2] = st
    q[0] = 0.0
    q[2] = ct

    for k in range(2, kmx + 1):
        if n < m:
            m = 0
            n = n + 1
            rr = rr * ratio
            fn = n
            gn = n - 1

        fm = m

        if m == n:
            if k != 3:
                one = math.sqrt(1 - 0.5 / fm)
                j = k - n - 1
                p[k - 1] = one * st * p[j - 1]
                q[k - 1] = one * (st * q[j - 1] + ct * p[j - 1])
                cl[m - 1] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0]
                sl[m - 1] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0]
        else:
            gmm = m ** 2
            one = math.sqrt(fn ** 2 - gmm)
            two = math.sqrt(gn ** 2 - gmm) / one
            three = (fn + gn) / one
            i = k - n
            j = i - n + 1
            p[k - 1] = three * ct * p[i - 1] - two * p[j - 1]
            q[k - 1] = three * (ct * q[i - 1] - st * p[i - 1]) - two * q[j - 1]

        lm = ll + l
        one = (tc * gh[lm - 1] + t * gh[lm + nc - 1]) * rr

        if m != 0:
            two = (tc * gh[lm] + t * gh[lm + nc]) * rr
            three = one * cl[m - 1] + two * sl[m - 1]
            x = x + three * q[k - 1]
            z = z - (fn + 1) * three * p[k - 1]

            if st != 0:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * fm * p[k - 1] / st
            else:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k - 1] * ct

            l = l + 2

        else:
            x = x + one * q[k - 1]
            z = z - (fn + 1) * one * p[k - 1]
            l = l + 1

        m = m + 1

    one = x
    x = x * Cd + z * sd
    z = z * Cd - one * sd

    # rotate to correct
    deg2rad = math.pi / 180
    phi = deg2rad * latitude_degrees
    lam = deg2rad * elongitude_degrees

    ecef_Q_ned = np.array(
        [
            [
                -math.sin(phi) * math.cos(lam),
                -math.sin(lam),
                -math.cos(phi) * math.cos(lam),
            ],
            [
                -math.sin(phi) * math.sin(lam),
                math.cos(lam),
                -math.cos(phi) * math.sin(lam),
            ],
            [math.cos(phi), 0.0, -math.sin(phi)],
        ]
    )
    print(np.dot(ecef_Q_ned, np.array([x, y, z])))


igrf13_5(date, latitude_degrees, elongitude_degrees, r_norm_km)
