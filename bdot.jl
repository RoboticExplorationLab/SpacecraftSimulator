using LinearAlgebra




function bdot_control_law(ᴺqᴮ,ω,max_dipoles,B_eci_T,eclipse)
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

    if eclipse
        m = zeros(3)
        return m
    else
        # attitude
        ᴺQᴮ = dcm_from_q(ᴺqᴮ)
        ᴮQᴺ = transpose(ᴺQᴮ)

        # magnetic field vector in body frame
        B_body_T = ᴮQᴺ*B_eci_T

        # bdot approximation
        bdot = -cross(ω,B_body_T)

        # bang-bang control law
        m = -max_dipoles .* sign.(bdot)

        return m
    end
end
