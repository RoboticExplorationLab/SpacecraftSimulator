using SatelliteToolbox, JLD2, LinearAlgebra

using JLD2
@load "/Users/kevintracy/devel/SpacecraftSimulator/rv_to_tle_SD/spacecraft_info.jld2" sc_info

tle_string = rv_to_tle([sc_info.jd_epc],[sc_info.r_eci],[sc_info.v_eci])
