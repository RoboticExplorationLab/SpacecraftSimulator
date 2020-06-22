using LinearAlgebra, SatelliteDynamics
const SD = SatelliteDynamics

# load in data types
include(joinpath(dirname(@__DIR__),"common/types.jl"))


"""Dynamics functions that are used in the propagator."""


function FODE(epc::Epoch,eci_state::Vec,u::Vec)::Vec

    # unpack state
    r_eci = eci_state[1:3]
    v_eci = eci_state[4:6]


    # spherical harmonic gravity stuff

    # ECI ECEF stuff
    # ECEF_Q_ECI = SD.rECItoECEF(epc)
    #
    # # acceleration
    # a_eci = SD.accel_gravity(r_eci,ECEF_Q_ECI,params.grav_deg,
    #                                           params.grav_order) + u

    # J2 only
    a_eci = FODE_J2(r_eci) + u


    return [v_eci;a_eci]
end

function FODE_J2(r_eci::Vec)::Vec

    # relavent parameters
    J2 = J2_EARTH
    μ = GM_EARTH

    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = (R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  (-μ/r^3)*[x*(1 - 1.5*J2*Re_r_sqr*(five_z_sqr - 1));
                      y*(1 - 1.5*J2*Re_r_sqr*(five_z_sqr - 1));
                      z*(1 - 1.5*J2*Re_r_sqr*(five_z_sqr - 3))]
end


function rk4_orbital(f::Function, t_n::RealorEpoch, x_n::Vec, u::Vec,
                                                              h::Real)::Vec
    """Runge-Kutta 4th order integration. Epoch for time.

    Args:
        ODE:           f(t,x,u)        | function
        initial time:  t_n             | epoch or float64
        initial cond:  x_n             | vec
        control input: u               | vec
        time step:     h               | scalar

    Returns:
        x_{n+1}
    """

    k1 = h*f(t_n,x_n,u)
    k2 = h*f(t_n+h/2,x_n+k1/2,u)
    k3 = h*f(t_n+h/2,x_n+k2/2,u)
    k4 = h*f(t_n+h,x_n+k3,u)

    # @infiltrate
    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end

function rk4_attitude(f::Function, t_n::RealorEpoch, x_n::Vec, m::Vec,
                                    B_eci_nT::Vec, τ::Vec, h::Real)::Vec
    """Runge-Kutta 4th order integration. Epoch for time.

    Args:
        ODE:                f(t,x,u)        | function
        initial time:       t_n             | epoch or float64
        initial cond:       x_n             | vec
        magnetic dipole:    m               | A⋅m²
        mag field vector:   B_eci_nT        | nT (eci)
        external torque:    τ               | N⋅m
        step size:          h               | s
    Returns:
        x_{n+1}
    """

    k1 = h*f(t_n,x_n,m,B_eci_nT,τ)
    k2 = h*f(t_n+h/2, x_n+k1/2, m, B_eci_nT, τ)
    k3 = h*f(t_n+h/2, x_n+k2/2, m, B_eci_nT, τ)
    k4 = h*f(t_n+h, x_n+k3, m, B_eci_nT, τ)

    x_np1 = (x_n + (1/6)*(k1+2*k2+2*k3 + k4));
    x_np1[1:4] = normalize(x_np1[1:4])
    # @infiltrate
    return x_np1
end

function spacecraft_eom(t, x, m, B_eci_nT, τ)
    """Spacecraft equations of motion.

    Args:
        t: time (epoch)
        x: state
            - quaternion (scalar last) ᴺqᴮ
            - angular velocity of the sc (rad/s, expressed in body)
        m: magnetic moment control (A⋅m^2)
        B_eci_nT: magnetic field in ECI (nT)
        τ: external torque (n⋅m^2)

    Returns:
        xdot: state derivative

    """
    # @infiltrate
    # error()
    # unpack state
    ᴺqᴮ = x[1:4]
    ᴺωᴮ = x[5:7]
    J = params.sc.J
    invJ = params.sc.invJ

    # magnetic field in the body frame
    ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))
    B_body_T = ᴮQᴺ*B_eci_nT

    # magnetic moment
    magnetic_moment = cross(m,B_body_T)

    # quaternion kinematics
    ᴺq̇ᴮ = 0.5 * ᴺqᴮ ⊙ [(ᴺωᴮ); 0.0]

    # angular acceleration
    ᴺαᴮ = invJ*(τ + magnetic_moment - ᴺωᴮ × (J * ᴺωᴮ) )

    # @infiltrate

    return [ᴺq̇ᴮ; ᴺαᴮ]
end

function EOE_eom(epc::Epoch, x::Vec, u::Vec,thruster_on::Bool)::Vec
    """Equations of motion for equinoctial elements in 2-body problem.

    Args:
        epc: time(::epoch)
        x: state
            - p: a*(1-e^2) semi-latus rectum       (m)
            - f: e*cos(ω + Ω)                      (non-dimensional)
            - g: e*sin(ω + Ω)                      (non-dimensional)
            - h: tan(i/2)*cos(Ω)                   (non-dimensional)
            - k: tan(i/2)*sin(Ω)                   (non-dimensional)
            - L: ω + Ω + θ argument of latitude    (rad)
            - mass: spacecraft mass                (kg)
        u: thruster acceleration                   (m/s^2)
        thruster_on: bool for thruster             ()

    Returns:
        ̇x: state derivative
    """

    # unpack state
    p,f,g,h,k,L,mass = x

    # throw error if the orbit goes hyperbolic
    if p<0
        error("orbit is hyperbolic")
    end

    # get ECI position and RTN_Q_ECI direction cosine matrix
    r_eci,v_eci = rv_from_eoe(x[1:6])
    RTN_Q_ECI = rtn_from_eci_mat(r_eci,v_eci)

    # equinoctial parameters for GVE's
    sinL = sin(L)
    cosL = cos(L)
    alpha_2 = h^2 - k^2
    s2 = 1+h^2 + k^2
    w = 1+f*cosL + g*sinL
    r = p/w

    # mass flow rate if thruster is on
    if thruster_on
        mass_dot = -params.thruster_mass_flow_rate
    else
        mass_dot = 0.0
    end

    # perturbation acceleration
    u_pert = zeros(3)

    # J2 acceleration only
    if params.J2_only
        u_pert += J2_from_EOE(p,f,g,h,k,L)
    end

    # spherical harmonic expansion gravity
    if params.spherical_harmonic_gravity
        u_pert += spherical_harmonic_gravity_model(epc, r_eci, RTN_Q_ECI)
    end

    # acceleration from drag (don't bother if altitude is > 800km)
    if params.drag && r < (params.Re + 800000)
        u_pert += drag_from_atmosphere(r_eci,v_eci,mass,epc, RTN_Q_ECI)
    end

    # sun 3rd body
    if params.sun_3rd_body
        u_pert += RTN_Q_ECI*accel_thirdbody_sun(epc,r_eci)
    end

    # moon 3rd body
    if params.moon_3rd_body
        u_pert += RTN_Q_ECI*accel_thirdbody_moon(epc,r_eci)
    end

    # add perturbing acceleration to input acceleration
    u_r, u_t, u_n = u_pert + u

    # precompute to save calculations
    sqrtpu = sqrt(p/params.mu)

    # GVE's for equinoctial orbital elements
    p_dot = (2*p/w)*sqrtpu*u_t
    f_dot = sqrtpu*( u_r*sinL    +    ((w+1)*cosL + f)*u_t/w    -
                                            (h*sinL - k*cosL)*g*u_n/w  )
    g_dot = sqrtpu*(-u_r*cosL   +    ((w+1)*sinL + g)*u_t/w    +
                                            (h*sinL - k*cosL)*f*u_n/w  )
    h_dot = sqrtpu*(s2*u_n/(2*w))*cosL
    k_dot = sqrtpu*(s2*u_n/(2*w))*sinL
    L_dot = sqrt(params.mu*p)*(w/p)^2 + (1/w)*sqrtpu*(h*sinL - k*cosL)*u_n

    # return state derivative
    return [p_dot;f_dot;g_dot;h_dot;k_dot;L_dot;mass_dot]
end

function spherical_harmonic_gravity_model(epc::Epoch, r_eci::Vec,
                                          RTN_Q_ECI::Mat)::Vec
    """Spherical harmonic gravity model up to 180 deg/order.

    Args:
        epc: epoch for time                     (Epoch)
        r_eci: cartesian position vector        (m)
        RTN_Q_ECI: DCM between RTN and ECI      ()
            [V_RTN = RTN_Q_ECI * V_ECI]

    Returns:
        a_rtn: acceleration in RTN frame        (m/s)

    Comments:
        The a_rtn acceleration is without the 2 body gravitational because
        the GVE's already account for it.
    """

    ECEF_Q_ECI = rECItoECEF(epc)

    a_eci = accel_gravity(r_eci,ECEF_Q_ECI,params.n_max,params.m_max)


    return a_eci
end

function J2_from_EOE(p::Real, f::Real, g::Real, h::Real, k::Real, L::Real)::Vec
    """J2 acceleration in RTN from equinoctial orbital elements.

    Args:
        p: a*(1-e^2) semi-latus rectum       (m)
        f: e*cos(ω + Ω)                      (non-dimensional)
        g: e*sin(ω + Ω)                      (non-dimensional)
        h: tan(i/2)*cos(Ω)                   (non-dimensional)
        k: tan(i/2)*sin(Ω)                   (non-dimensional)
        L: ω + Ω + θ argument of latitude    (rad)

    Returns:
        a_rtn: acceleration in the RTN frame  (m/s)
    """

    J2 =  params.J2

    # constants for EOE
    mu = params.mu
    w = 1+f*cos(L) + g*sin(L)
    r = p/w
    Re = params.Re

    # precompute trig functions
    sinL = sin(L)
    cosL = cos(L)

    # accelerations
    a_r = -((3*mu*J2*Re^2)/(2*r^4))*(1 -   (  (12*(h*sinL-k*cosL)^2 )
                                                    /   ((1+h^2 + k^2)^2))  )
    a_t = -((12*mu*J2*Re^2)/(r^4))*((  ((h*sinL - k*cosL)*(h*cosL+k*sinL) )
                                                    /  ((1+h^2 + k^2)^2))  )
    a_n = -((6*mu*J2*Re^2)/(r^4))*((  ((h*sinL - k*cosL)*(1-h^2-k^2) )
                                                    /  ((1+h^2 + k^2)^2))  )
    return [a_r;a_t;a_n]
end

function drag_from_atmosphere(r_eci::Vec, v_eci::Vec, mass::Real, epc::Epoch,
                              RTN_Q_ECI::Mat)::Vec
    """Atmospheric drag acceleration in RTN.

    Args:
        r_eci: position vector in ECI         (m)
        v_eci: velocity vector in ECI         (m/s)
        mass: spacecraft mass                 (kg)
        epc: spacecraft time                  (Epoch)
        RTN_Q_ECI: DCM between rtn and eci    ()

    Returns:
        drag acceleration in rtn              (m/s^2)
    """

    # ECEF_Q_ECI DCM
    ECEF_Q_ECI = rECItoECEF(epc)

    # s/c position in ecef
    r_ecef = ECEF_Q_ECI*r_eci

    # geodetic coordinates
    geod = sECEFtoGEOD(r_ecef)

    # atmosphere model
    rho = density_nrlmsise00(epc,geod)

    # relative velocity
    v_rel = v_eci - transpose(ECEF_Q_ECI)*cross(params.w_earth,r_ecef)

    # acceleration from relative velocity and air density
    a_drag_eci = -.5*rho*params.S_area*params.C_d*v_rel*norm(v_rel)/mass

    return a_drag_eci
end


function eclipse_and_srp(epc::Epoch, r_eci::Vec, mass::Real,
                         RTN_Q_ECI::Mat)::Tuple{Bool,Vec}
    """Eclipse check and SRP acceleration.

    Args:
        epc: time                            (Epoch)
        r_eci: position vector eci           (m)
        mass: s/c mass                       (kg)
        RTN_Q_ECI: DCM between basis         ()

    Returns:
        eclipse: bool, true if eclipse       ()
        a_srp: srp acceleration in RTN       (m/s^2)
    """

    # get the sun position
    r_sun_eci = sun_position(epc)

    # normalize sun position
    e_sun_eci = normalize(r_sun_eci)

    #amount of r_eci in the direction of the sun
    proj = dot(r_eci,e_sun_eci)

    # if this is positive, we aren't in eclipse
    if proj > 0
        eclipse =  false
    else
        # part of r_eci orthogonal to sun vector
        ortho_r_eci = r_eci - proj*e_sun_eci

        # if the spacecraft is outside the shadow of earth
        if norm(ortho_r_eci)> params.Re
            eclipse = false
        else
            eclipse = true
        end
    end
    return eclipse
end



function eclipse_and_srp(epc::Epoch, r_eci::Vec, mass::Real,
                         RTN_Q_ECI::Mat)::Tuple{Bool,Vec}
    """Eclipse check and SRP acceleration.

    Args:
        epc: time                            (Epoch)
        r_eci: position vector eci           (m)
        mass: s/c mass                       (kg)
        RTN_Q_ECI: DCM between basis         ()

    Returns:
        eclipse: bool, true if eclipse       ()
        a_srp: srp acceleration in RTN       (m/s^2)
    """

    # get the sun position
    r_sun_eci = sun_position(epc)

    # normalize sun position
    e_sun_eci = normalize(r_sun_eci)

    #amount of r_eci in the direction of the sun
    proj = dot(r_eci,e_sun_eci)

    # if this is positive, we aren't in eclipse
    if proj > 0
        eclipse =  false
    else
        # part of r_eci orthogonal to sun vector
        ortho_r_eci = r_eci - proj*e_sun_eci

        # if the spacecraft is outside the shadow of earth
        if norm(ortho_r_eci)> params.Re
            eclipse = false
        else
            eclipse = true
        end
    end

    # if params.srp is set to true, and the spacecraft is not in eclipse
    if params.srp && !eclipse
        # SRP force vector
        # normalize(r_eci - r_sun_eci) is the unit vector from sun to s/c
        F_srp = params.P_sun*params.cr*params.A_srp*(normalize(r_eci - r_sun_eci))

        # SRP acceleration vector in rtn
        return eclipse, RTN_Q_ECI*F_srp/mass
    else
        # otherwise return 0 for a_srp_rtn
        return eclipse, zeros(3)
    end

end

function eclipse_and_srp(epc::Epoch, r_eci::Vec, mass::Real,
                         RTN_Q_ECI::Mat)::Tuple{Bool,Vec}
    """Eclipse check and SRP acceleration.

    Args:
        epc: time                            (Epoch)
        r_eci: position vector eci           (m)
        mass: s/c mass                       (kg)
        RTN_Q_ECI: DCM between basis         ()

    Returns:
        eclipse: bool, true if eclipse       ()
        a_srp: srp acceleration in RTN       (m/s^2)
    """

    # get the sun position
    r_sun_eci = sun_position(epc)

    # normalize sun position
    e_sun_eci = normalize(r_sun_eci)

    #amount of r_eci in the direction of the sun
    proj = dot(r_eci,e_sun_eci)

    # if this is positive, we aren't in eclipse
    if proj > 0
        eclipse =  false
    else
        # part of r_eci orthogonal to sun vector
        ortho_r_eci = r_eci - proj*e_sun_eci

        # if the spacecraft is outside the shadow of earth
        if norm(ortho_r_eci)> params.Re
            eclipse = false
        else
            eclipse = true
        end
    end

    # if params.srp is set to true, and the spacecraft is not in eclipse
    if params.srp && !eclipse
        # SRP force vector
        # normalize(r_eci - r_sun_eci) is the unit vector from sun to s/c
        F_srp = params.P_sun*params.cr*params.A_srp*(normalize(r_eci - r_sun_eci))

        # SRP acceleration vector in rtn
        return eclipse, RTN_Q_ECI*F_srp/mass
    else
        # otherwise return 0 for a_srp_rtn
        return eclipse, zeros(3)
    end

end
