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

    xtraj = [np.zeros(Nx) for _ in range(N)]
    xtraj[0] = x0
    utraj = [np.zeros(Nu) for _ in range(N - 1)]

    # # Forward simulate using initial controls
    # for k in range(0, N - 1):
    #     # J = (
    #     #     J
    #     #     + 0.5 * (xtraj[:, k] - xg) @ Q @ (xtraj[:, k] - xg)
    #     #     + 0.5 * R * utraj[:, k] @ utraj[:, k]
    #     # )
    #     J += cost(Q, R, (xtraj[:, k] - xg).transpose(), utraj[:, k].transpose())
    #     (xtraj[:, k + 1]) = rkstep_xdot(xtraj[:, k], utraj[:, k], dt)

    # # J = J + 0.5 * (xtraj[:, N - 1] - xg) @ Qf @ (xtraj[:, N - 1] - xg)
    # J += 0.5 * quad(Qf, (xtraj[:, N - 1] - xg).transpose())

    # print("their way")
    # print(J)
    # print("easy way")
    J = (N - 1) * 0.5 * quad(Q, x0 - xg) + 0.5 * quad(Qf, x0 - xg)
    # pdb.set_trace()
    Jhist.append(J)

    S = np.zeros((Nx, Nx))
    s = np.zeros(Nx)
    K = np.zeros((Nu, Nx * (N - 1)))
    # l = (tol + 1) * np.ones((Nu, N - 1))
    l = [np.zeros(Nu) for _ in range(N - 1)]

    count = 0
    for ii in range(100):

        count += 1

        S = Qf
        # s = Qf @ (xtraj[:, N - 1] - xg)
        s = mul2(Qf, (xtraj[N - 1] - xg).transpose())

        # Backward pass
        for k in range(N - 2, -1, -1):

            # Calculate cost gradients for this time step
            # q = Q @ (xtraj[:, k] - xg)
            # r = R * utraj[:, k]
            q = mul2(Q, (xtraj[k] - xg).transpose())
            r = mul2(R, utraj[k].transpose())
            # pdb.set_trace()
            # Make assignments for ease of reading
            # Ak = A[:, Nx * k : Nx * (k + 1)]
            # Bk = B[:, Nu * k : Nu * (k + 1)]
            Ak, Bk = rkstep_jacobians(xtraj[k], utraj[k], dt)

            # Calculate l and K
            # LH = R + Bk.T @ S @ Bk
            invLH = np.linalg.inv(R + mul3(Bk.transpose(), S, Bk))

            # LH = R + quad(S, Bk)
            l[k] = mul2(invLH, (r + mul2(Bk.transpose(), s)))
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
                + mul3(Kk.transpose(), R, l[k])
                + mul2((Ak - mul2(Bk, Kk)).transpose(), (s - mul3(S, Bk, l[k])))
            )
            S = Snew
            s = snew

        # Forward pass line search with new l and K
        # unew = np.zeros((Nu, N - 1))
        # xnew = np.zeros((Nx, N))
        # xnew[:, 0] = x0
        xnew = [np.zeros(Nx) for _ in range(N)]
        xnew[0] = x0
        unew = [np.zeros(Nu) for _ in range(N - 1)]
        alpha = 1.0
        Jnew = J + 1
        while Jnew > J:
            Jnew = 0
            for k in range(0, N - 1):
                unew[k] = (
                    utraj[k]
                    - alpha * l[k]
                    - mul2(K[:, Nx * k : Nx * (k + 1)], (xnew[k] - xtraj[k]))
                )
                (xnew[k + 1]) = rkstep_xdot(xnew[k], unew[k], dt)

                Jnew += cost(Q, R, (xnew[k] - xg).transpose(), unew[k].transpose())

            # Jnew = Jnew + 0.5 * (xnew[:, N - 1] - xg).T @ Qf @ (xnew[:, N - 1] - xg)
            Jnew += 0.5 * quad(Qf, (xnew[N - 1] - xg).transpose())
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
        print(2 * alpha)
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


# def rkstep_xdot(x0, u0, dt):

#     # Define constants
#     Nx = x0.shape[0]
#     Nu = u0.shape[0]

#     xdot1 = dynamics_xdot(0, x0, u0)
#     xdot2 = dynamics_xdot(0, x0 + 0.5 * xdot1 * dt, u0)

#     x1 = x0 + dt * xdot2

#     return x1


def rkstep_jacobians(x1, u0, dt):

    # # Define constants
    # Nx = x0.shape[0]
    # Nu = u0.shape[0]

    # A1, B1 = dynamics_jacobians(0, x0, u0)
    # xdot1 = dynamics_xdot(0, x0, u0)
    # A2, B2 = dynamics_jacobians(0, x0 + 0.5 * xdot1 * dt, u0)

    # A = np.eye(2) + dt * A2 + 0.5 * dt * dt * A2 @ A1
    # B = dt * B2 + 0.5 * dt * dt * A2 @ B1

    # return A, B
    Nx = 2

    # x1 = x0
    t = 0.0

    xdot1, A1, B1 = pendulumDynamics_py(t, x1, u0)
    k1 = xdot1 * dt

    x2 = x1 + 0.5 * k1
    xdot2, A2, B2 = pendulumDynamics_py(t + dt * 0.5, x2, u0)
    k2 = dt * xdot2

    x3 = x1 + 0.5 * k2
    xdot3, A3, B3 = pendulumDynamics_py(t + dt * 0.5, x3, u0)
    k3 = dt * xdot3

    x4 = x1 + k3
    xdot4, A4, B4 = pendulumDynamics_py(t + dt, x4, u0)
    # k4 = dt * xdot4

    # x_tp1 = x1 + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    # RK4 method
    # A_d
    dk1_dx1 = dt * A1
    dx2_dx1 = np.eye(Nx) + 0.5 * dk1_dx1
    dk2_dx1 = dt * mul2(A2, dx2_dx1)
    dx3_dx1 = np.eye(Nx) + 0.5 * dk2_dx1
    dk3_dx1 = dt * mul2(A3, dx3_dx1)
    dx4_dx1 = np.eye(Nx) + dk3_dx1
    dk4_dx1 = dt * mul2(A4, dx4_dx1)
    A_d = np.eye(Nx) + (1 / 6) * (dk1_dx1 + 2 * dk2_dx1 + 2 * dk3_dx1 + dk4_dx1)

    # B_d
    dk1_du = dt * B1
    dx2_du = 0.5 * dk1_du
    dk2_du = dt * mul2(A2, dx2_du) + dt * B2
    dx3_du = 0.5 * dk2_du
    dk3_du = dt * mul2(A3, dx3_du) + dt * B3
    dx4_du = 1.0 * dk3_du
    dk4_du = dt * mul2(A4, dx4_du) + dt * B4
    B_d = (1 / 6) * (dk1_du + 2 * dk2_du + 2 * dk3_du + dk4_du)

    return A_d, B_d


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


def rkstep_xdot(x1, u0, dt):

    # xdot1 = dynamics_xdot(0, x0, u0)
    # k1 = xdot1 * dt
    # x2
    # xdot2 = dynamics_xdot(0, x0 + 0.5 * xdot1 * dt, u0)

    # x1 = x0 + dt * xdot2
    t = 0.0

    xdot1 = dynamics_xdot(t, x1, u0)
    k1 = xdot1 * dt

    x2 = x1 + 0.5 * k1
    xdot2 = dynamics_xdot(t + dt * 0.5, x2, u0)
    k2 = dt * xdot2

    x3 = x1 + 0.5 * k2
    xdot3 = dynamics_xdot(t + dt * 0.5, x3, u0)
    k3 = dt * xdot3

    x4 = x1 + k3
    xdot4 = dynamics_xdot(t + dt, x4, u0)
    k4 = dt * xdot4

    return x1 + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


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

    # xtraj = np.concatenate(xtraj, axis=0)
    # numpy.stack( LIST, axis=0 )
    # pdb.set_trace()
    xtraj = np.stack(xtraj, axis=0).T
    utraj = np.stack(utraj, axis=0).T
    # Plot results
    fig, ax = plt.subplots(2, 2)
    fig.suptitle("iLQR results")
    ax[0, 0].plot(xtraj[0, :])
    ax[0, 0].set_title("angle")
    ax[0, 1].plot(xtraj[1, :])
    ax[0, 1].set_title("angular rate")
    ax[1, 0].plot(utraj[0, :])
    ax[1, 0].set_title("Control")
    ax[1, 1].plot(Jhist)
    ax[1, 1].set_title("Cost")
    plt.show()


if __name__ == "__main__":
    main()

