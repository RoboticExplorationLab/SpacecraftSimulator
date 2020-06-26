function attitude_loop!(orbital_state,attitude_state,kk,B_eci,r_sun_eci,
                        eclipse,inner_loop_t_vec,attitude_state_sense,
                        B_body_sense_vec,s_body_sense_vec,epc_orbital,dt_attitude)


# attitude dynamics inner loop
for jj = 1:length(inner_loop_t_vec)-1

    # index of current step n
    index_n = (kk-1)*(length(inner_loop_t_vec)-1) + jj

    # true state
    ᴺqᴮ_true = attitude_state[1:4,index_n]
    ω_true = attitude_state[5:7,index_n]
    r_eci_true = orbital_state[1:3,kk]
    v_eci_true = orbital_state[4:6,kk]

    # environmental stuff
    B_eci_T_true = B_eci[:,kk]
    r_sun_eci_true = r_sun_eci
    eclipse_true = eclipse

    # ---------------truth measurements-----------------
    # attitude
    ᴺQᴮ_true = dcm_from_q(ᴺqᴮ_true)
    ᴮQᴺ_true = transpose(ᴺQᴮ_true)

    # magnetic field vector in the body
    B_body_T_true = ᴮQᴺ_true*B_eci_T_true

    # sun flux
    I_vec_true = sun_flux(r_sun_eci_true,r_eci_true,ᴺQᴮ_true)

    # sun vector in body
    s_body_true = s_body_from_I(I_vec_true)

    # eclipse sense
    eclipse_sense = eclipse_true

    # sense and actuate opposite steps

    #---------------sensor measurements-----------------
    ᴺqᴮ_sense = ᴺqᴮ_true
    ᴺQᴮ_sense = dcm_from_q(ᴺqᴮ_sense)
    ᴮQᴺ_sense = transpose(ᴺQᴮ_sense)
    ω_sense = ω_true
    B_body_T_sense = B_body_T_true
    s_body_sense = s_body_true

    # store it all in the main arrays
    attitude_state_sense[1:4,index_n] = ᴺqᴮ_sense
    attitude_state_sense[5:7,index_n] = ω_sense
    B_body_sense_vec[1:3,index_n] = B_body_T_sense
    s_body_sense_vec[1:3,index_n] = s_body_sense

    if isodd(index_n)
        # store it all in the main arrays
        attitude_state_sense[1:4,index_n] = ᴺqᴮ_sense
        attitude_state_sense[5:7,index_n] = ω_sense
        B_body_sense_vec[1:3,index_n] = B_body_T_sense
        s_body_sense_vec[1:3,index_n] = s_body_sense

        # control
        sc_mag_moment = zeros(3)
    else
        attitude_state_sense[1:7,index_n] = attitude_state_sense[1:7,(index_n-1)]
        B_body_sense_vec[1:3,index_n] = B_body_sense_vec[1:3,(index_n-1)]
        s_body_sense_vec[1:3,index_n] = s_body_sense_vec[1:3,(index_n-1)]

        ω_sense = attitude_state_sense[5:7,index_n]
        B_body_T_sense = B_body_sense_vec[1:3,index_n]

        #---------------control law-------------------------
        # bdot Control law (julia)
        sc_mag_moment = bdot_control_law(ω_sense,
                                         params.sc.max_dipoles,
                                         B_body_T_sense,
                                         eclipse_sense)

    end

    # disturbance torques
    τ = zeros(3)

    # update the attitude with RK4
    attitude_state[:,index_n+1] =rk4_attitude(spacecraft_eom,
    epc_orbital, attitude_state[:,index_n], sc_mag_moment, B_eci_T_true, τ, dt_attitude)

end

end
