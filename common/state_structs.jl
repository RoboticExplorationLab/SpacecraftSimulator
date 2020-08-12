"""Functions relating to the state struct."""

mutable struct truth_state_struct
    r_eci           :: Array{Array{Float64,1},1}
    v_eci           :: Array{Array{Float64,1},1}
    ᴺqᴮ             :: Array{Array{Float64,1},1}
    ᴺQᴮ             :: Array{Array{Float64,2},1} # derived
    ω               :: Array{Array{Float64,1},1}
    B_eci           :: Array{Array{Float64,1},1} # derived
    eclipse_hist    :: Array{Bool,1}             # derived
    r_sun_eci       :: Array{Array{Float64,1},1} # derived
    I_sun_flux      :: Array{Array{Float64,1},1} # derived
    sun_body        :: Array{Array{Float64,1},1} # derived
    B_body          :: Array{Array{Float64,1},1} # derived
    epc_orbital     :: Epoch
end

function orbital_truth_struct_update!(truth::truth_state_struct,k::Int)
    """Fill in the truth struct with derived properties for orbital changes."""

    # sun position
    truth.r_sun_eci[k] = SD.sun_position(truth.epc_orbital)

    # eclipse
    truth.eclipse_hist[k] = eclipse_check(truth.r_eci[k],truth.r_sun_eci[k])

    # ECI magnetic field (T)
    truth.B_eci[k] = IGRF13(truth.r_eci[k],truth.epc_orbital)

end

function attitude_truth_struct_update!(truth::truth_state_struct,k::Int,jj::Int)
    """Fill in the truth struct with derived properties for attitude changes."""

    # pull attitude and get the DCM
    truth.ᴺQᴮ[jj]   = dcm_from_q(truth.ᴺqᴮ[jj])

    # magnetic field in the body
    truth.B_body[jj] = transpose(truth.ᴺQᴮ[jj])*truth.B_eci[k]

    # Sun Flux
    truth.I_sun_flux[jj] = sun_flux(truth.r_sun_eci[k],truth.r_eci[k],
                                   truth.ᴺQᴮ[jj],truth.eclipse_hist[k])

    # sun_body
    truth.sun_body[jj] = sun_body_normalized(truth.r_sun_eci[k],truth.r_eci[k],
                                            truth.ᴺQᴮ[jj])

end


function initialize_struct(struct_type_name   ::DataType,
                           time_params        ::NamedTuple,
                           initial_conditions ::NamedTuple)

    """Initialize the truth struct with pre-allocated arrays and IC's."""
    t_vec_orbital  = time_params.t_vec_orbital
    t_vec_attitude = time_params.t_vec_attitude

    orbital3 = fill(zeros(3),length(t_vec_orbital))
    orbital6 = fill(zeros(6),length(t_vec_orbital))
    attitude3 = fill(zeros(3),length(t_vec_attitude))
    attitude4 = fill(zeros(4),length(t_vec_attitude))
    attitude3x3 = fill(zeros(3,3),length(t_vec_attitude))
    attitude6 = fill(zeros(6),length(t_vec_attitude))
    attitude7 = fill(zeros(7),length(t_vec_attitude))
    orbital_bool = fill(false,length(t_vec_orbital))
    epc_orbital = Epoch("2018-12-20")

    truth = struct_type_name(copy(orbital3),      # r_eci
                             copy(orbital3),      # v_eci
                             copy(attitude4),     # q
                             copy(attitude3x3),   # Q
                             copy(attitude3),     # ω
                             copy(orbital3),      # B_eci
                             copy(orbital_bool),  # eclipse
                             copy(orbital3),      # r_sun_eci
                             copy(attitude6),     # I_sun_flux
                             copy(attitude3),     # sun_body
                             copy(attitude3),     # B_body
                             epc_orbital          # time
                             )

    # initial conditions
    truth.r_eci[1]  = initial_conditions.eci_rv_0[1:3]
    truth.v_eci[1]  = initial_conditions.eci_rv_0[4:6]
    truth.ᴺqᴮ[1] =  initial_conditions.ᴺqᴮ0
    truth.ω[1] = initial_conditions.ω0
    truth.epc_orbital = initial_conditions.epc_orbital
        return truth
end



function mat_from_vec(a::Array{Array{Float64,1},1})

    rows = length(a[1])
    columns = length(a)
    A = zeros(rows,columns)

    for i = 1:columns
        A[:,i] = a[i]
    end

    return A
end
