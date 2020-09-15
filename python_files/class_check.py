import numpy as np
import math


class trajoptclass:
    gh_igrf13 = [
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

    def __init__(self, N, dt, Nx, Nu):
        self.K = [np.zeros((Nu, Nx)) for _ in range(N - 1)]
        self.l = [np.zeros(Nu) for _ in range(N - 1)]
        self.xtraj = [np.zeros(Nx) for _ in range(N)]
        self.utraj = [np.zeros(Nu) for _ in range(N - 1)]
        self.xnew = [np.zeros(Nx) for _ in range(N)]
        self.unew = [np.zeros(Nu) for _ in range(N - 1)]

    def testfx(self):
        print(self.xnew)
        print("hey")
        print(self.gh_igrf13)


TrajOpt = trajoptclass(25, 12.0, 6, 3)

TrajOpt.testfx()

print(TrajOpt.gh_igrf13)
