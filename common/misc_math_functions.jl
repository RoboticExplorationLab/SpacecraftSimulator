using LinearAlgebra
using Random
# load in data types
include(joinpath(dirname(@__DIR__),"common/types.jl"))

"""Miscellaneous math functions"""


function create_ellipse(height::Number, width::Number, xc::Number,
            yc::Number)::Tuple{Vec,Vec}
    """Creates x and y coordinates for an ellipse of given size.

    Args:
        height: total height :: Number
        width:  total width :: Number
        xc:     x coordinate of ellipse center :: Number
        yc:     y coordinate of ellipse center :: Number

    Returns:
        x_hist: vector of x's :: AbstractArray{Float64,1}
        y_hist: vector of y's :: AbstractArray{Float64,1}
    """

    # get vec of angles
    theta_vec = 0:.01:2 * pi

    # open up arrays for x and y
    x_hist = zeros(length(theta_vec))
    y_hist = zeros(length(theta_vec))

    # populate with values (offset π/2 so it starts at [0,-ymax])
    for i = 1:length(theta_vec)
        x_hist[i] = cos(theta_vec[i] - pi / 2)
        y_hist[i] = sin(theta_vec[i] - pi / 2)
    end

    # scale to height and width, and move to have center at (xc,yc)
    x_hist = x_hist * width / 2 .+ xc
    y_hist = y_hist * height / 2 .+ yc

    return x_hist, y_hist
end


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

    return vec([mat[3, 2]; mat[1, 3]; mat[2, 1]])
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
