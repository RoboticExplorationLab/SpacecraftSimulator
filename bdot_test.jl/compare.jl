using LinearAlgebra




function bdot_control_law(ω,max_dipoles,B_body_T,eclipse)
    """Bdot control law for detumbling.

    Args:
        ᴺqᴮ: attitude quaternion
        ω: angular velocity of the spacecraft wrt ECI, expressed in the body
        max_dipoles: max magnetic dipole from spacecraft, (A⋅m²)
        B_body_T: Earth magnetic field vector expressed in the body (T)

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

        # bdot approximation
        bdot = -cross(ω,B_body_T)

        # bang-bang control law
        m = -max_dipoles .* tanh.(1e4*bdot)

        return m
    end
end

ω = [.2,-.4,-.1]
max_dipoles = [8.8e-3,1.373e-2,8.2e-3]
B_body_T = [4e-5,.5e-5,9.4e-5]
eclipse = false

@show bdot_control_law(ω,max_dipoles,B_body_T,eclipse)
