using SatelliteDynamics
"""Commonly used data types."""

# 1d array that must be real
Vec = AbstractArray{<:Real,1}

# vector of vectors
VecofVecs = Union{ Array{Array{Float64,1},1}, Array{Array{Float32,1},1} }

# vector of matrices
VecofMats = Union{Array{Array{Float64,2},1},Array{Array{Float32,2},1}}

# 2d array that must be real
Mat = AbstractArray{<:Real,2}

# either a number or an Epoch (SatelliteDynamics.jl)
RealorEpoch = Union{Epoch,Real}
