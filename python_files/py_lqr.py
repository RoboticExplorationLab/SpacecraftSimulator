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


def iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol):
    """
	Simple iLQR for testing C++ functions
	"""
    Nx = x0.shape[0]
    Nu = utraj0.shape[0]
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
        J = (
            J
            + 0.5 * (xtraj[:, k] - xg) @ Q @ (xtraj[:, k] - xg)
            + 0.5 * R * utraj[:, k] @ utraj[:, k]
        )
        (
            xtraj[:, k + 1],
            A[:, Nx * k : Nx * (k + 1)],
            B[:, Nx * k : Nx * (k + 1)],
        ) = rkstep_py(xtraj[:, k], utraj[:, k], dt)

    J = J + 0.5 * (xtraj[:, N - 1] - xg) @ Qf @ (xtraj[:, N - 1] - xg)
    Jhist.append(J)

    S = np.zeros((Nx, Nx))
    s = np.zeros(Nx)
    K = np.zeros((Nu, Nx * (N - 1)))
    l = (tol + 1) * np.ones((Nu, N - 1))

    count = 0
    while np.amax(np.abs(l)) > tol:

        count += 1

        S = Qf
        s = Qf @ (xtraj[:, N - 1] - xg)

        # Backward pass
        for k in range(N - 2, -1, -1):

            # Calculate cost gradients for this time step
            q = Q @ (xtraj[:, k] - xg)
            r = R * utraj[:, k]

            # Make assignments for ease of reading
            Ak = A[:, Nx * k : Nx * (k + 1)]
            Bk = B[:, Nu * k : Nu * (k + 1)]

            # Calculate l and K
            LH = R + Bk.T @ S @ Bk
            l[:, k] = np.linalg.solve(LH, (r + Bk.T @ s))
            K[:, Nx * k : Nx * (k + 1)] = np.linalg.solve(LH, Bk.T @ S @ Ak)

            # Calculate new S and s
            Kk = K[:, Nx * k : Nx * (k + 1)]
            Snew = Q + R * Kk.T @ Kk + (Ak - Bk @ Kk).T @ S @ (Ak - Bk @ Kk)
            snew = (
                q
                - Kk.T @ r
                + R * Kk.T @ l[:, k]
                + (Ak - Bk @ Kk).T @ (s - S @ Bk @ l[:, k])
            )
            S = Snew
            s = snew

        # Forward pass line search with new l and K
        unew = np.zeros((Nu, N - 1))
        xnew = np.zeros((Nx, N))
        xnew[:, 0] = xtraj[:, 0]
        alpha = 1.0
        Jnew = J + 1
        while Jnew > J:
            Jnew = 0
            for k in range(0, N - 1):
                unew[:, k] = (
                    utraj[:, k]
                    - alpha * l[:, k]
                    - K[:, Nx * k : Nx * (k + 1)] @ (xnew[:, k] - xtraj[:, k])
                )
                (
                    xnew[:, k + 1],
                    A[:, Nx * k : Nx * (k + 1)],
                    B[:, Nx * k : Nx * (k + 1)],
                ) = rkstep_py(xnew[:, k], unew[:, k], dt)
                Jnew = (
                    Jnew
                    + 0.5 * (xnew[:, k] - xg).T @ Q @ (xnew[:, k] - xg)
                    + 0.5 * R * unew[:, k].T @ unew[:, k]
                )

            Jnew = Jnew + 0.5 * (xnew[:, N - 1] - xg).T @ Qf @ (xnew[:, N - 1] - xg)
            alpha = 0.5 * alpha

        xtraj = xnew
        utraj = unew
        J = Jnew
        Jhist.append(J)

        print("Iteration {}".format(count))
        print("Final l = {}".format(np.max(np.abs(l))), "alpha = {}".format(2 * alpha))

    return xtraj, utraj, K, Jhist


def rkstep_py(x0, u0, dt):

    # Define constants
    Nx = x0.shape[0]
    Nu = u0.shape[0]

    xdot1, A1, B1 = pendulumDynamics_py(0, x0, u0)
    xdot2, A2, B2 = pendulumDynamics_py(0, x0 + 0.5 * xdot1 * dt, u0)

    x1 = x0 + dt * xdot2

    # A1 = dxdot1[:, 0:Nx]
    # A2 = dxdot2[:, 0:Nx]
    # B1 = dxdot1[:, Nx:]
    # B2 = dxdot2[:, Nx:]

    A = np.eye(2) + dt * A2 + 0.5 * dt * dt * A2 @ A1
    B = dt * B2 + 0.5 * dt * dt * A2 @ B1

    return x1, A, B


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


# Test the algorithm
def main():

    N = 750
    Nx = 2
    Nu = 1

    Qf = np.eye(Nx) * 30
    Q = np.eye(Nx) * 0.01
    R = 0.03

    x0 = np.zeros(2)
    xg = np.array([np.pi, 0])
    utraj0 = np.zeros((Nu, N - 1))

    dt = 0.01
    tol = 0.2

    xtraj, utraj, K, Jhist = iLQRsimple_py(x0, xg, utraj0, Q, R, Qf, dt, tol)

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

