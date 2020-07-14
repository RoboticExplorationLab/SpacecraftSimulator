using LinearAlgebra
using Random
# load in data types
# include(joinpath(dirname(@__DIR__),"common/types.jl"))

"""Miscellaneous math functions"""


function wrap_to_2pi(theta::Real)::Real
    """Takes an angle theta and returns the same angle theta ∈ [0,2π].

    Args:
        theta: angle in radians :: Float64

    Returns:
        theta: angle in radians wrapped :: Float64
    """

    # if angle is negative
    if theta < 0.0
        theta = -2*pi*(abs(theta)/(2*pi)-floor(abs(theta/(2*pi)))) + 2*pi

    # if angle is positive
    else
        theta = 2*pi*(abs(theta)/(2*pi)-floor(abs(theta/(2*pi))))
    end

    return theta
end

function wrap_to_pm_pi(theta::Real)::Real
    """Takes an angle theta and returns the same angle theta ∈ [-π,π].

    Args:
        theta: angle in radians :: Float64

    Returns:
        theta: angle in radians wrapped :: Float64

    Comments:
        'pm' in the function title stands for 'plus minus', or ±
    """

    # wrap the angle to 2π
    theta = wrap_to_2pi(theta)

    # if the angle is > pi, subtract 2π s.t. it's between ±π
    if theta>pi
        theta -= 2*pi
    end

    return theta
end

function skew_from_vec(v::Vec)::Mat
    """Skew symmetric matrix from a vector.

    Summary:
        Takes 3 component vector and returns matrix after skew_from_vec
        operator has been applied. This is the same as the cross product
        matrix. cross(a,b) = skew_from_vec(a)*b

    Args:
        v::Vector

    Returns:
        Skew symmetric matrix from the given vector :: AbstractArray{Float64,2}
    """

    v = float(vec(v))
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end

function vec_from_skew(mat::Mat)::Vec
    """Converts 3x3 skew symmetric matrix to a 3x1 vector.

    Args:
        mat: 3x3 skew symmetric matrix

    Returns:
        3x1 vector
    """

    return [mat[3, 2]; mat[1, 3]; mat[2, 1]]
end

function dcm_from_phi(phi::Vec)::Mat
    """DCM from axis angle (phi)"""
    return skew_expm(skew_from_vec(phi))
end


function skew_expm(B::Mat)::Mat
    """Expm for skew symmetric matrices.

    Summary:
        matrix exponential for skew symmetric matrices. about 40% faster than
        the standard exp.jl function. This function can be used to take the
        skew symmetric skew_from_vec matrix of an axis angle vector, and
        producing the orthogonal Direction Cosine Matrix (DCM)
        skew_from_vec corresponds.

    Args:
        B: skew symmetric matrix :: AbstractArray{Float64,2}

    Returns:
        orthognal matrix :: AbstractArray{Float64,2}
    """

    # axis (from axis-angle)
    phi = vec_from_skew(B)

    # angle (from axis-angle)
    theta = norm(phi)

    # axis
    if theta == 0
        r = [0.0; 0.0; 0.0]
    else
        r = phi / theta
    end

    # closed form skew symmetric matrix exponential
    return (I + sin(theta) * skew_from_vec(r) + (1.0 - cos(theta)) *
    skew_from_vec(r) * skew_from_vec(r))
end

function ortho_logm(Q::Mat)::Mat
    """Matrix logarithm for 3x3 orthogonal matrices (like DCM's).

    Summary:
        This is both faster and more robust than log.jl (180 degree rotations)

    Args:
        Q: orthogonal matrix (like a DCM) :: AbstractArray{Float64,2}

    Returns:
        skew symmetric matrix :: AbstractArray{Float64,2}
    """

    val = (tr(Q) - 1) / 2

    if abs(val - 1) < 1e-10
        # no rotation
        phi = [0.0; 0.0; 0.0]
    elseif abs(val + 1) < 1e-10
        # 180 degree rotation
        M = I + Q
        r = M[1, :] / norm(M[1, :])
        theta = pi
        phi = r * theta
    else
        # normal rotation (0-180 degrees)
        theta = acos(val)
        r = -(1 / (2 * sin(theta))) *
            [Q[2, 3] - Q[3, 2]; Q[3, 1] - Q[1, 3]; Q[1, 2] - Q[2, 1]]
        phi = r * theta
    end

    return skew_from_vec(phi)
end

function rand_in_range(lower::Real, upper::Real)::Real
    """Random number within range with a uniform distribution.

    Args:
        lower: lower value in range
        upper: upper value in range

    Returns:
        random number in the specified range
    """

    delta = upper - lower
    return lower + rand() * delta
end


function active_rotation_axis_angle(axis::Vec,theta::Real,vector::Vec)::Vec
    """Actively rotate a vector using an axis and angle.

    Args:
        axis: axis of rotation (unit norm)
        theta: angle to rotate (rad)

    Returns:
        rotated vector (unit norm)
    """

    # check if axis is unit norm
    if !isapprox(norm(axis),1.0,rtol=1e-6)
        error("axis is not unit norm")
    end

    return skew_expm(skew_from_vec(theta*axis))*vector
end


function ⊙(q1, q2)
    """Quaternion multiplication, hamilton product, scalar last"""

    v1 = q1[1:3]
    s1 = q1[4]
    v2 = q2[1:3]
    s2 = q2[4]

    return [(s1 * v2 + s2 * v1 + cross(v1, v2));(s1 * s2 - dot(v1, v2))]

end

function dcm_from_q(q)
    """DCM from quaternion, hamilton product, scalar last"""

    # pull our the parameters from the quaternion
    q1,q2,q3,q4 = q

    # DCM
    Q = [(2*q1^2+2*q4^2-1)   2*(q1*q2 - q3*q4)   2*(q1*q3 + q2*q4);
          2*(q1*q2 + q3*q4)  (2*q2^2+2*q4^2-1)   2*(q2*q3 - q1*q4);
          2*(q1*q3 - q2*q4)   2*(q2*q3 + q1*q4)  (2*q3^2+2*q4^2-1)]

    return Q
end

function qconj(q)
    """Conjugate of the quaternion (scalar last)"""

    return [-q[1:3]; q[4]]
end

function phi_from_q(q)
    """axis angle from quaternion (scalar last)"""

    v = q[1:3]
    s = q[4]
    normv = norm(v)

    if normv == 0.0
        return zeros(3)
    else
        r = v / normv
        θ = 2 * atan(normv, s)
        return r * θ
    end
end

function q_from_phi(ϕ)
    """Quaternion from axis angle vector, scalar last"""

    θ = norm(ϕ)
    if abs(θ) < 0.0000000001
        return [0; 0; 0; 1.0]
    else
        r = ϕ / θ
        return [r * sin(θ / 2); cos(θ / 2)]
    end
end

function q_from_dcm(dcm)
    """Kane/Levinson convention, scalar last"""
    R = dcm
    T = R[1,1] + R[2,2] + R[3,3]
    if T> R[1,1] && T > R[2,2] && T>R[3,3]
        q4 = .5*sqrt(1+T)
        r  = .25/q4
        q1 = (R[3,2] - R[2,3])*r
        q2 = (R[1,3] - R[3,1])*r
        q3 = (R[2,1] - R[1,2])*r
    elseif R[1,1]>R[2,2] && R[1,1]>R[3,3]
        q1 = .5*sqrt(1-T + 2*R[1,1])
        r  = .25/q1
        q4 = (R[3,2] - R[2,3])*r
        q2 = (R[1,2] + R[2,1])*r
        q3 = (R[1,3] + R[3,1])*r
    elseif R[2,2]>R[3,3]
        q2 = .5*sqrt(1-T + 2*R[2,2])
        r  = .25/q2
        q4 = (R[1,3] - R[3,1])*r
        q1 = (R[1,2] + R[2,1])*r
        q3 = (R[2,3] + R[3,2])*r
    else
        q3 = .5*sqrt(1-T + 2*R[3,3])
        r  = .25/q3
        q4 = (R[2,1] - R[1,2])*r
        q1 = (R[1,3] + R[3,1])*r
        q2 = (R[2,3] + R[3,2])*r
    end
    q = [q1;q2;q3;q4]
    if q4<0
        q = -q
    end

    return q
end

function randq()
    return normalize(randn(4))
end

function q_shorter(q)

    if q[4]<0
        q = -q
    end
    return q
end


function H()
    """matrix for converting vector to pure quaternion. Scalar last"""
    return [1 0 0;
            0 1 0;
            0 0 1;
            0 0 0.0]
end
