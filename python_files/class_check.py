import numpy as np
import math


class trajopt:
    def __init__(self, N, dt, Nx, Nu):
        self.K = [np.zeros((Nu, Nx)) for _ in range(N - 1)]
        self.l = [np.zeros(Nu) for _ in range(N - 1)]
        self.xtraj = [np.zeros(Nx) for _ in range(N)]
        self.utraj = [np.zeros(Nu) for _ in range(N - 1)]
        self.xnew = [np.zeros(Nx) for _ in range(N)]
        self.unew = [np.zeros(Nu) for _ in range(N - 1)]

    def testfx(self):
        print(self.xnew)


TrajOpt = trajopt(25, 12.0, 6, 3)

TrajOpt.testfx()
