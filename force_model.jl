




function force_model(eci_state::Vec,epc_orbital::Epoch,r_sun_eci::Vec)::Vec
    """Orbital force modeling"""
    
    r_eci =  eci_state[1:3]
    v_eci =  eci_state[4:6]

    # sun position and eclipse check
    # r_sun_eci = SD.sun_position(epc_orbital)
    # eclipse = eclipse_check(r_eci, r_sun_eci)

    a_perturbation = zeros(3)

    # atmospheric drag
    ρ = density_harris_priester(r_eci, r_sun_eci)
    ecef_Q_eci = SD.rECItoECEF(epc_orbital)
    a_drag = accel_drag(orbital_state[:,kk], ρ, params.sc.mass,
                        params.sc.area, params.sc.cd, ecef_Q_eci)
    a_perturbation += a_drag

    # third body accel from sun
    a_3rd_body_sun = accel_thirdbody_sun(r_eci,r_sun_eci)
    a_perturbation += a_3rd_body_sun

    # third body accel from moon
    a_3rd_body_moon = accel_thirdbody_moon(epc_orbital,r_eci)
    a_perturbation += a_3rd_body_moon

    # accel from srp
    a_srp = accel_srp(r_eci,r_sun_eci,params.sc.mass,params.sc.area)
    a_perturbation += a_srp

    return a_perturbation
end
