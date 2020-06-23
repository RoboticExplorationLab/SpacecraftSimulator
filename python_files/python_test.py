import numpy as np


def foo(inn):

    return np.array([1 * inn, 2 * inn, 3 * inn])


def bdot_control_law(q, w, max_dipoles, B_eci_T, eclipse):
    """Bdot control law for detumbling.

    Args:
        ᴺqᴮ: attitude quaternion
        ω: angular velocity of the spacecraft wrt ECI, expressed in the body
        max_dipoles: max magnetic dipole from spacecraft, (A⋅m²)
        B_eci_T: Earth magnetic field vector expressed in ECI (T)

    Returns:
        m: spacecraft magnetic dipole (A⋅m²)

    Ref:
        Fundamentals of Spacecraft Attitude Determination and Control (7.5.1)
        F. Landis Markley, John L. Crassidis
    """

    if eclipse:
        m = np.zeros(3)
        return m
    else:
        # attitude
        n_Q_b = dcm_from_q(q)
        b_Q_n = n_Q_b.T

        # magnetic field vector in body frame
        B_body_T = b_Q_n @ B_eci_T

        # bdot approximation
        bdot = -np.cross(w, B_body_T)

        # bang-bang control law
        m = -max_dipoles * np.sign(bdot)

        return m


def dcm_from_q(q):
    """DCM from quaternion, hamilton product, scalar last"""

    # pull our the parameters from the quaternion
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    q4 = q[3]

    # DCM
    Q = np.array(
        [
            [
                (2 * q1 ** 2 + 2 * q4 ** 2 - 1),
                2 * (q1 * q2 - q3 * q4),
                2 * (q1 * q3 + q2 * q4),
            ],
            [
                2 * (q1 * q2 + q3 * q4),
                (2 * q2 ** 2 + 2 * q4 ** 2 - 1),
                2 * (q2 * q3 - q1 * q4),
            ],
            [
                2 * (q1 * q3 - q2 * q4),
                2 * (q2 * q3 + q1 * q4),
                (2 * q3 ** 2 + 2 * q4 ** 2 - 1),
            ],
        ]
    )

    return Q
