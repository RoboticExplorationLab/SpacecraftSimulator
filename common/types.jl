using SatelliteDynamics
"""Commonly used data types."""

# 1d array that must be real
Vec = AbstractArray{<:Real,1}

# 2d array that must be real
Mat = AbstractArray{<:Real,2}

# either a number or an Epoch (SatelliteDynamics.jl)
RealorEpoch = Union{Epoch,Real}
