from __future__ import division, print_function, absolute_import

import numpy as np


def fAndG(A, b, x):
    assert isinstance(A, np.ndarray)
    dim = A.shape
    assert len(dim) == 2
    A_rows = dim[0]
    A_cols = dim[1]
    assert isinstance(b, np.ndarray)
    dim = b.shape
    assert len(dim) == 1
    b_rows = dim[0]
    assert isinstance(x, np.ndarray)
    dim = x.shape
    assert len(dim) == 1
    x_rows = dim[0]
    assert b_rows == A_rows
    assert x_rows == A_cols

    t_0 = (A).dot(x) - b
    functionValue = np.sum(np.log(t_0))
    gradient = (A.T).dot((np.ones(b_rows) / t_0))

    return functionValue, gradient


def checkGradient(A, b, x):
    # numerical gradient checking
    # f(x + t * delta) - f(x - t * delta) / (2t)
    # should be roughly equal to inner product <g, delta>
    t = 1e-6
    delta = np.random.randn(3)
    f1, _ = fAndG(A, b, x + t * delta)
    f2, _ = fAndG(A, b, x - t * delta)
    f, g = fAndG(A, b, x)
    print(
        "approximation error",
        np.linalg.norm((f1 - f2) / (2 * t) - np.tensordot(g, delta, axes=1)),
    )


def generateRandomData():
    A = np.random.randn(3, 3)
    b = np.random.randn(3)
    x = np.random.randn(3)

    return A, b, x


if __name__ == "__main__":
    A, b, x = generateRandomData()
    functionValue, gradient = fAndG(A, b, x)
    print("functionValue = ", functionValue)
    print("gradient = ", gradient)

    print("numerical gradient checking ...")
    checkGradient(A, b, x)
