import numpy as np
import math


def mul2(a, b):
    return np.dot(a, b)


def mul3(a, b, c):
    return np.dot(mul2(a, b), c)


def quad(Q, x):
    return mul3(x.transpose(), Q, x)


Q = np.array(
    [
        [0.32166826, 0.08742232, 0.13397362],
        [0.152078, 0.84912493, 0.97621681],
        [0.56134641, 0.14855114, 0.61592457],
    ]
)

x = np.array([0.45313514, 0.70570864, 0.36871091])

print(quad(Q, x.transpose()))

