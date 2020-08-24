import numpy as np


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
        -21.1,
        -11.7,  # 2020
        15.3,
        14.9,
        13.7,
        3.6,
        -16.5,
        -6.9,  # 2020
        -0.3,
        2.8,
        5.0,
        8.4,
        -23.4,
        2.9,  # 2020
        11.0,
        -1.5,
        9.8,
        -1.1,
        -5.1,
        -13.2,  # 2020
        -6.3,
        1.1,
        7.8,
        8.8,
        0.4,
        -9.3,
    ]
)


def colatd_from_latd(latd):
    return 90 - latd


date = 2020.3

R_EARTH = 6.378136300e6  # m

lat = np.rad2deg(33)
elong = np.rad2deg(76)
alt = 1.05 * R_EARTH / 1000
order = 6

colat = colatd_from_latd(lat)


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

# gh = gh[3256:end]
# ll  = 3255
ll = 0
nmx = order

nc = int(nmx * (nmx + 2))
# nc = 195

kmx = int((nmx + 1) * (nmx + 2) / 2)
# kmx = 105

# allocate
cl = np.zeros(nmx)
sl = np.zeros(nmx)
p = np.zeros(kmx)
q = np.zeros(kmx)


r = alt
ct = np.cos(colat * np.pi / 180)
st = np.sin(colat * np.pi / 180)
cl[0] = np.cos(elong * np.pi / 180)
sl[0] = np.sin(elong * np.pi / 180)
Cd = 1.0
sd = 0.0
l = 1
m = 1
n = 0

ratio = 6371.2 / r
rr = ratio ** 2

p[1 - 1] = 1.0
p[3 - 1] = st
q[1 - 1] = 0.0
q[3 - 1] = ct

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
            one = np.sqrt(1 - 0.5 / fm)
            j = k - n - 1
            p[k - 1] = one * st * p[j - 1]
            q[k - 1] = one * (st * q[j - 1] + ct * p[j - 1])
            cl[m - 1] = cl[m - 1 - 1] * cl[1 - 1] - sl[m - 1 - 1] * sl[1 - 1]
            sl[m - 1] = sl[m - 1 - 1] * cl[1 - 1] + cl[m - 1 - 1] * sl[1 - 1]
    else:
        gmm = m ** 2
        one = np.sqrt(fn ** 2 - gmm)
        two = np.sqrt(gn ** 2 - gmm) / one
        three = (fn + gn) / one
        i = k - n
        j = i - n + 1
        p[k - 1] = three * ct * p[i - 1] - two * p[j - 1]
        q[k - 1] = three * (ct * q[i - 1] - st * p[i - 1]) - two * q[j - 1]

    lm = ll + l
    one = (tc * gh[lm - 1] + t * gh[lm + nc - 1]) * rr

    if m != 0:
        two = (tc * gh[lm + 1 - 1] + t * gh[lm + nc + 1 - 1]) * rr
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

print([x, y, z])

