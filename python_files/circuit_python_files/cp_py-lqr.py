import time
import board
import math
import ulab as np


def mul2(a, b):
    """Multiply two arrays"""
    return np.linalg.dot(a, b)


def mul3(a, b, c):
    """Multiply three arrays"""
    return np.linalg.dot(mul2(a, b), c)


def mul4(a, b, c, d):
    """Multiply four arrays"""
    return np.linalg.dot(mul3(a, b, c), d)


def mul5(a, b, c, d, e):
    """Multiply five arrays"""
    return np.linalg.dot(mul4(a, b, c, d), e)


def quad(Q, x):
    """Quadratic x'*Q*x"""
    return mul3(x.transpose(), Q, x)


def cost(Q, R, x, u):
    """Quadratic cost function .5*x'Qx + .5*u'Ru"""
    return (0.5 * quad(Q, x) + 0.5 * quad(R, u))[0]


def iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol):
    """Simple iLQR with no initial guess"""
    Nx = 2
    Nu = 1
    N = utraj0.shape[1] + 1

    J = 0.0
    Jhist = []

    xtraj = [np.zeros(Nx) for _ in range(N)]
    xtraj[0] = x0
    utraj = [np.zeros(Nu) for _ in range(N - 1)]

    # initial cost of an all zeros guess
    J = (N - 1) * 0.5 * quad(Q, (x0 - xg).transpose())[0] + 0.5 * quad(
        Qf, (x0 - xg).transpose()
    )[0]

    Jhist.append(J)
    S = np.zeros((Nx, Nx))
    s = np.zeros(Nx)
    K = np.zeros((Nu, Nx * (N - 1)))
    K = [np.zeros((Nu, Nx)) for _ in range(N - 1)]
    l = [np.zeros(Nu) for _ in range(N - 1)]

    count = 0
    for ii in range(50):

        count += 1

        S = Qf
        s = mul2(Qf, (xtraj[N - 1] - xg).transpose())

        # Backward pass
        for k in range(N - 2, -1, -1):

            # Calculate cost gradients for this time step
            q = mul2(Q, (xtraj[k] - xg).transpose())
            r = mul2(R, utraj[k].transpose())

            Ak, Bk = rkstep_jacobians(xtraj[k], utraj[k], dt)

            # Calculate l and K
            invLH = np.linalg.inv(R + mul3(Bk.transpose(), S, Bk))

            l[k] = mul2(invLH, (r + mul2(Bk.transpose(), s)))
            K[k] = mul2(invLH, mul3(Bk.transpose(), S, Ak))

            # Calculate new S and s

            Snew = Q + mul3(K[k].transpose(), R, K[k]) + quad(S, (Ak - mul2(Bk, K[k])))

            snew = (
                q
                - mul2(K[k].transpose(), r)
                + mul3(K[k].transpose(), R, l[k])
                + mul2((Ak - mul2(Bk, K[k])).transpose(), (s - mul3(S, Bk, l[k])))
            )
            S = Snew
            s = snew

        # Forward pass line search with new l and K
        xnew = [np.zeros(Nx) for _ in range(N)]
        xnew[0] = x0
        unew = [np.zeros(Nu) for _ in range(N - 1)]
        alpha = 1.0
        Jnew = J + 1
        while Jnew > J:
            Jnew = 0.0
            for k in range(0, N - 1):
                unew[k] = (
                    utraj[k]
                    - alpha * l[k]
                    - mul2(K[k], (xnew[k] - xtraj[k]).transpose())
                )
                (xnew[k + 1]) = rkstep_xdot(xnew[k], unew[k], dt)

                Jnew += cost(Q, R, (xnew[k] - xg).transpose(), unew[k].transpose())

            Jnew += 0.5 * quad(Qf, (xnew[N - 1] - xg).transpose())[0]
            alpha = 0.5 * alpha

        dJ = J - Jnew
        xtraj = xnew
        utraj = unew
        J = Jnew
        Jhist.append(J)

        if dJ < 0.1:
            break

        print(ii)
        print(2 * alpha)
        print(J)

    return xtraj, utraj, K, Jhist


def rkstep_xdot(x0, u0, dt):

    xdot1 = dynamics_xdot(0, x0, u0)
    xdot2 = dynamics_xdot(0, x0 + 0.5 * xdot1 * dt, u0)

    x1 = x0 + dt * xdot2

    return x1


def rkstep_jacobians(x0, u0, dt):

    # Define constants

    A1, B1 = dynamics_jacobians(0, x0, u0)
    xdot1 = dynamics_xdot(0, x0, u0)
    A2, B2 = dynamics_jacobians(0, x0 + 0.5 * xdot1 * dt, u0)

    A = np.eye(2) + dt * A2 + 0.5 * dt * dt * mul2(A2, A1)
    B = dt * B2 + 0.5 * dt * dt * mul2(A2, B1)

    return A, B


def dynamics_xdot(t, x, u):
    Nx = 2

    m = 1.0
    b = 0.1
    lc = 0.5
    I = 0.25
    g = 9.81

    xdot = np.zeros(Nx)

    xdot[0] = x[1]
    xdot[1] = (u - m * g * lc * math.sin(x[0]) - b * x[1]) / I

    return xdot


def dynamics_jacobians(t, x, u):
    Nx = 2

    m = 1.0
    b = 0.1
    lc = 0.5
    I = 0.25
    g = 9.81

    A = np.array([[0, 1], [-m * g * lc * math.cos(x[0]) / I, -b / I]])
    B = np.array([[0], [1 / I]])

    return A, B


# Test the algorithm
def main2():

    N = 75
    Nx = 2
    Nu = 1

    Qf = np.eye(Nx) * 30
    Q = np.eye(Nx) * 0.01
    R = np.eye(1) * 0.03

    x0 = np.zeros(2)
    xg = np.array([math.pi, 0])
    utraj0 = np.zeros((Nu, N - 1))

    dt = 0.1
    tol = 0.35

    xtraj, utraj, K, Jhist = iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol)


main2()
