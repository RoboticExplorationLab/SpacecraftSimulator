"""
 * Author: Nick Goodson 
 * Jan 18th 2020
 *
 * This code is a template, different sections can be pulled out and replaced
 * with calls to C++ functions to test the C++ functionality
 *
"""

import pdb
import numpy as np
import matplotlib.pyplot as plt

# filename = inspect.getframeinfo(inspect.currentframe()).filename
# currentdir = os.path.dirname(os.path.abspath(filename))
# parentdir = os.path.dirname(currentdir)
# gncdir = os.path.dirname(parentdir)

# sys.path.insert(0, parentdir)
# sys.path.insert(0, gncdir)

# import CPP modules here
# import iLQRsimple_cpp as ilqr


def mul2(a, b):
    return np.dot(a, b)


def mul3(a, b, c):
    return np.dot(mul2(a, b), c)


def mul4(a, b, c, d):
    return np.dot(mul3(a, b, c), d)


def mul5(a, b, c, d, e):
    return np.dot(mul4(a, b, c, d), e)


def quad(Q, x):
    return mul3(x.transpose(), Q, x)


def cost(Q, R, x, u):
    return 0.5 * quad(Q, x) + 0.5 * quad(R, u)


def iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol):
    """
	Simple iLQR for testing C++ functions
	"""
    Nx = 2
    Nu = 1
    N = utraj0.shape[1] + 1

    A = np.zeros((Nx, Nx * (N - 1)))
    B = np.zeros((Nx, Nu * (N - 1)))

    J = 0
    Jhist = []
    xtraj = np.zeros((Nx, N))
    xtraj[:, 0] = x0
    utraj = utraj0

    # Forward simulate using initial controls
    for k in range(0, N - 1):
        # J = (
        #     J
        #     + 0.5 * (xtraj[:, k] - xg) @ Q @ (xtraj[:, k] - xg)
        #     + 0.5 * R * utraj[:, k] @ utraj[:, k]
        # )
        J += cost(Q, R, (xtraj[:, k] - xg).transpose(), utraj[:, k].transpose())
        (xtraj[:, k + 1]) = rkstep_xdot(xtraj[:, k], utraj[:, k], dt)

    # J = J + 0.5 * (xtraj[:, N - 1] - xg) @ Qf @ (xtraj[:, N - 1] - xg)
    J += 0.5 * quad(Qf, (xtraj[:, N - 1] - xg).transpose())

    # pdb.set_trace()
    Jhist.append(J)

    S = np.zeros((Nx, Nx))
    s = np.zeros(Nx)
    K = np.zeros((Nu, Nx * (N - 1)))
    l = (tol + 1) * np.ones((Nu, N - 1))

    count = 0
    for ii in range(100):

        count += 1

        S = Qf
        # s = Qf @ (xtraj[:, N - 1] - xg)
        s = mul2(Qf, (xtraj[:, N - 1] - xg).transpose())

        # Backward pass
        for k in range(N - 2, -1, -1):

            # Calculate cost gradients for this time step
            # q = Q @ (xtraj[:, k] - xg)
            # r = R * utraj[:, k]
            q = mul2(Q, (xtraj[:, k] - xg).transpose())
            r = mul2(R, utraj[:, k].transpose())
            # pdb.set_trace()
            # Make assignments for ease of reading
            # Ak = A[:, Nx * k : Nx * (k + 1)]
            # Bk = B[:, Nu * k : Nu * (k + 1)]
            Ak, Bk = rkstep_jacobians(xtraj[:, k], utraj[:, k], dt)

            # Calculate l and K
            # LH = R + Bk.T @ S @ Bk
            invLH = np.linalg.inv(R + mul3(Bk.transpose(), S, Bk))

            # LH = R + quad(S, Bk)
            l[:, k] = mul2(invLH, (r + mul2(Bk.transpose(), s)))
            K[:, Nx * k : Nx * (k + 1)] = mul2(invLH, mul3(Bk.transpose(), S, Ak))

            # Calculate new S and s
            Kk = K[:, Nx * k : Nx * (k + 1)]
            # Snew = Q + R * Kk.T @ Kk + (Ak - Bk @ Kk).T @ S @ (Ak - Bk @ Kk)
            # pdb.set_trace()
            Snew = Q + mul3(Kk.transpose(), R, Kk) + quad(S, (Ak - mul2(Bk, Kk)))
            # snew = (
            #     q
            #     - Kk.T @ r
            #     + R * Kk.T @ l[:, k]
            #     + (Ak - Bk @ Kk).T @ (s - S @ Bk @ l[:, k])
            # )
            snew = (
                q
                - mul2(Kk.transpose(), r)
                + mul3(Kk.transpose(), R, l[:, k])
                + mul2((Ak - mul2(Bk, Kk)).transpose(), (s - mul3(S, Bk, l[:, k])))
            )
            S = Snew
            s = snew

        # Forward pass line search with new l and K
        unew = np.zeros((Nu, N - 1))
        xnew = np.zeros((Nx, N))
        xnew[:, 0] = x0
        alpha = 1.0
        Jnew = J + 1
        while Jnew > J:
            Jnew = 0
            for k in range(0, N - 1):
                unew[:, k] = (
                    utraj[:, k]
                    - alpha * l[:, k]
                    - mul2(K[:, Nx * k : Nx * (k + 1)], (xnew[:, k] - xtraj[:, k]))
                )
                (xnew[:, k + 1]) = rkstep_xdot(xnew[:, k], unew[:, k], dt)

                Jnew += cost(
                    Q, R, (xnew[:, k] - xg).transpose(), unew[:, k].transpose()
                )

            # Jnew = Jnew + 0.5 * (xnew[:, N - 1] - xg).T @ Qf @ (xnew[:, N - 1] - xg)
            Jnew += 0.5 * quad(Qf, (xnew[:, N - 1] - xg).transpose())
            alpha = 0.5 * alpha

        dJ = J - Jnew
        xtraj = xnew
        utraj = unew
        J = Jnew
        Jhist.append(J)

        if dJ < 0.1:
            break

        # print("Iteration {}".format(count))
        # print("Final l = {}".format(np.max(np.abs(l))), "alpha = {}".format(2 * alpha))
        print(ii)
        print(alpha)
        print(J)

    return xtraj, utraj, K, Jhist


# def rkstep_py(x0, u0, dt):

#     # Define constants
#     Nx = x0.shape[0]
#     Nu = u0.shape[0]

#     xdot1, A1, B1 = pendulumDynamics_py(0, x0, u0)
#     xdot2, A2, B2 = pendulumDynamics_py(0, x0 + 0.5 * xdot1 * dt, u0)

#     x1 = x0 + dt * xdot2

#     # A1 = dxdot1[:, 0:Nx]
#     # A2 = dxdot2[:, 0:Nx]
#     # B1 = dxdot1[:, Nx:]
#     # B2 = dxdot2[:, Nx:]

#     A = np.eye(2) + dt * A2 + 0.5 * dt * dt * A2 @ A1
#     B = dt * B2 + 0.5 * dt * dt * A2 @ B1

#     return x1, A, B


def rkstep_xdot(x0, u0, dt):

    # Define constants
    Nx = x0.shape[0]
    Nu = u0.shape[0]

    xdot1 = dynamics_xdot(0, x0, u0)
    xdot2 = dynamics_xdot(0, x0 + 0.5 * xdot1 * dt, u0)

    x1 = x0 + dt * xdot2

    return x1


def rkstep_jacobians(x0, u0, dt):

    # Define constants
    Nx = x0.shape[0]
    Nu = u0.shape[0]

    A1, B1 = dynamics_jacobians(0, x0, u0)
    xdot1 = dynamics_xdot(0, x0, u0)
    A2, B2 = dynamics_jacobians(0, x0 + 0.5 * xdot1 * dt, u0)

    A = np.eye(2) + dt * A2 + 0.5 * dt * dt * A2 @ A1
    B = dt * B2 + 0.5 * dt * dt * A2 @ B1

    return A, B


def pendulumDynamics_py(t, x, u):
    Nx = x.shape[0]

    m = 1.0
    b = 0.1
    lc = 0.5
    I = 0.25
    g = 9.81

    xdot = np.zeros(Nx)
    xdot[0] = x[1]
    xdot[1] = (u - m * g * lc * np.sin(x[0]) - b * x[1]) / I

    A = np.array([[0, 1], [-m * g * lc * np.cos(x[0]) / I, -b / I]])
    B = np.array([[0], [1 / I]])
    # dxdot = np.hstack((A, B))

    return xdot, A, B


def dynamics_xdot(t, x, u):
    Nx = x.shape[0]

    m = 1.0
    b = 0.1
    lc = 0.5
    I = 0.25
    g = 9.81

    xdot = np.zeros(Nx)
    xdot[0] = x[1]
    xdot[1] = (u - m * g * lc * np.sin(x[0]) - b * x[1]) / I

    return xdot


def dynamics_jacobians(t, x, u):
    Nx = x.shape[0]

    m = 1.0
    b = 0.1
    lc = 0.5
    I = 0.25
    g = 9.81

    A = np.array([[0, 1], [-m * g * lc * np.cos(x[0]) / I, -b / I]])
    B = np.array([[0], [1 / I]])
    # dxdot = np.hstack((A, B))

    return A, B


# Test the algorithm
def main():

    N = 75
    Nx = 2
    Nu = 1

    Qf = np.eye(Nx) * 30
    Q = np.eye(Nx) * 0.01
    R = np.eye(1) * 0.03

    x0 = np.zeros(2)
    xg = np.array([np.pi, 0])
    utraj0 = np.zeros((Nu, N - 1))

    dt = 0.1
    tol = 0.35

    xtraj, utraj, K, Jhist = iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol)

    # # Plot results
    # fig, ax = plt.subplots(2, 2)
    # fig.suptitle("iLQR results")
    # ax[0, 0].plot(xtraj[0, :])
    # ax[0, 0].set_title("angle")
    # ax[0, 1].plot(xtraj[1, :])
    # ax[0, 1].set_title("angular rate")
    # ax[1, 0].plot(utraj[0, :])
    # ax[1, 0].set_title("Control")
    # ax[1, 1].plot(Jhist)
    # ax[1, 1].set_title("Cost")
    # plt.show()


if __name__ == "__main__":
    main()

