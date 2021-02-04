using LinearAlgebra



function clamp_bdot!(bdot)
    for i = 1:3
        if bdot[i]>0
            bdot[i]=min(.5*bdot[i],1)
        else
            bdot[i]=max(.5*bdot[i],-1)
        end
    end
end


# function test_clamp_bdot()
#     x = -10:.1:10
#     y = zeros(length(x))
#     for i = 1:length(x)
#         bdot = x[i]*ones(3)
#         clamp_bdot!(bdot)
#         y[i] = bdot[1]
#     end
#
#     mat"
#     figure
#     hold on
#     plot($x,$y)
#     hold off
#     "
# end
function bdot_control_law(ω,max_dipoles,B_body_T,eclipse)
    """Bdot control law for detumbling.

    Args:
        ᴺqᴮ: attitude quaternion
        ω: angular velocity of the spacecraft wrt ECI, expressed in body (rad)
        max_dipoles: max magnetic dipole from spacecraft, (A⋅m²)
        B_body_T: Earth magnetic field vector expressed in the body (T)

    Returns:
        m: spacecraft magnetic dipole (A⋅m²)

    Ref:
        Fundamentals of Spacecraft Attitude Determination and Control (7.5.1)
        F. Landis Markley, John L. Crassidis
    """


    # bdot approximation (NOTE: scaled it up)
    bdot = -1e5*cross(ω,B_body_T)

    clamp_bdot!(bdot)

    # tanh control law
    m = (!eclipse)*(-max_dipoles .* bdot)

    return m

end

function testit()
ω = [.8,-.4,-.1]
max_dipoles = [8.8e-3,1.373e-2,8.2e-3]
B_body_T = [4e-5,.5e-5,9.4e-5]
eclipse = false

m = bdot_control_law(ω,max_dipoles,B_body_T,eclipse)
@show m
end
testit()
