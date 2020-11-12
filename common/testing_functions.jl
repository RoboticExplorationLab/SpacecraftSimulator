"""Functions to be used in unit testing."""
# load in data types
# include(joinpath(dirname(@__DIR__),"common/types.jl"))

function isbetween(val::Float64,lower::Real,upper::Real) ::Bool
    """Testing function to determining if a value is with a specified range.

    Args:
        val: value to be tested
        lower: lower value in range
        upper: upper value in range

    returns:
        boolean that is true if val is within range [lower,upper]
    """

    if val>=lower && val<=upper
        return true
    else
        return false
    end
end

function matrix_isapprox(A::Mat, B::Mat,atol::Real)::Bool
    """Testing function for isapprox for two matrices.

    Args:
        A: a matrix
        B: a matrix

    Returns:
        boolean that is true if matrices A ≈ B
    """

    if size(A) != size(B)
        return false
    else
        if maximum(abs.(A-B)) > atol
            return false
        else
            return true
        end
    end
end


function angle_isapprox(a::Real,b::Real,tol::Real)::Bool
    """Unit testing function for two angles.

    Args:
        a: first angle    (rad)
        b: second angle   (rad)

    Returns:
        isapprox bool: true if ∠a ≈ ∠b, false if not

    Comments:
        This function is used for when the angles represent the same physical
        angle but have different values. An example would be θ₁ = 0.0 and
        θ₂ = 2π-.00001. The values are different but the angle is close.
    """

    # test the sin and cos of each angle
    test_bool_sin = isapprox(sin(a),sin(b),rtol = tol)
    test_bool_cos = isapprox(cos(a),cos(b),rtol = tol)

    # if both are equal then a and b represent the same angle
    if test_bool_sin && test_bool_cos
        return true
    else
        return false
    end
end
