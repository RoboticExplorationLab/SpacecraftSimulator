using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox
using YAML
const SD = SatelliteDynamics
const ST = SatelliteToolbox

function config(path_to_config_yaml)
      # load in the config.yml
      # data = YAML.load(open(joinpath(dirname(@__DIR__),"sim/config.yml")))
      data = YAML.load(open(joinpath(dirname(dirname(@__FILE__)),path_to_config_yaml)))

      # break it into the sections
      sc_properties = data["sc_properties"]
      gravity = data["gravity"]
      environment = data["environment"]
      init_oe = data["initial_orbital_elements"]
      propagator_parameters = data["propagator_parameters"]
      initial_attitude_conditions = data["initial_attitude_conditions"]

      # spacecraft inertia (expressed in body axes)
      # J = [1.959e-4 2016.333e-9 269.176e-9;
      #     2016.333e-9 1.999e-4 2318.659e-9;
      #     269.176e-9 2318.659e-9 1.064e-4]

      J = [sc_properties["inertia_rows"][1]';
           sc_properties["inertia_rows"][2]';
           sc_properties["inertia_rows"][3]']
      # J = [6400 -76.4 -25.6;-76.4 4730 -40;-25.6 -40 8160]
      # J = diagm([1;2;3.0])
      invJ = inv(J)

      # spacecraft mass
      # sc_mass = 177.808e-3
      # sc_area = .01
      # sc_cd = 1.5
      faces = [sc_properties["faces"][1]';
               sc_properties["faces"][2]';
               sc_properties["faces"][3]';
               sc_properties["faces"][4]';
               sc_properties["faces"][5]';
               sc_properties["faces"][6]']
      # face = vcat(faces[1]',faces[2]',faces[3]',)

      sc_mass = sc_properties["mass"]
      sc_area = sc_properties["area"]
      sc_cd = sc_properties["cd"]

      # max_dipoles = [8.8e-3;1.373e-2;8.2e-3] #amps

      max_dipoles = sc_properties["max_dipoles"]

      sc = (mass = sc_mass, area = sc_area, cd = sc_cd,J=J,invJ=invJ,
            max_dipoles = max_dipoles,faces=faces)

      # spherical harmonic expansion
      spherical_harmonic_gravity_bool = gravity["spherical_harmonic_gravity_bool"]
      grav_deg = gravity["gravity_degree"]
      grav_order = gravity["gravity_order"]

      # timing stuff
      epc_orbital = Epoch(propagator_parameters["epoch"])
      dt_orbital = propagator_parameters["dt_orbital"]
      dt_attitude = propagator_parameters["dt_attitude"]
      dt_controller = propagator_parameters["dt_controller"]
      tf = propagator_parameters["final_time_days"]*24*3600     # seconds

      # time vectors for simulation
      t_vec_orbital = 0:dt_orbital:tf
      inner_loop_t_vec = 0:dt_attitude:dt_orbital
      t_vec_attitude = 0:dt_attitude:(t_vec_orbital[end]+dt_orbital-dt_attitude)

      # time parameters named tuple
      time_params = (dt_orbital = dt_orbital, dt_attitude = dt_attitude,
                     dt_controller = dt_controller,tf = tf,
                     t_vec_orbital=t_vec_orbital, inner_loop_t_vec=inner_loop_t_vec,
                     t_vec_attitude = t_vec_attitude)

      # initial conditions
      ᴺqᴮ0 = initial_attitude_conditions["q0"]
      ω0 = initial_attitude_conditions["w0"]


      # orbital elements
      sma = init_oe["sma"]
      ecc = init_oe["ecc"]
      inc = init_oe["inc"]
      raan = init_oe["raan"]
      argp = init_oe["argp"]
      M = init_oe["M"]

      oe0 = [sma;ecc;inc;raan;argp;M]

      # convert to ECI r and v
      eci_rv_0 = sOSCtoCART(oe0, use_degrees=false)

      initial_conditions = (epc_orbital = epc_orbital,ᴺqᴮ0=ᴺqᴮ0,ω0=ω0,
                            eci_rv_0 = eci_rv_0,oe0 = oe0)


      params = (sc=sc, time_params=time_params,initial_conditions=initial_conditions,
                grav_deg = grav_deg, grav_order = grav_order,
                spherical_harmonic_gravity_bool = spherical_harmonic_gravity_bool)

      return params,initial_conditions, time_params
end



# read the TLE's
# line1 = "1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991"
# line2 = "2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102"
#
#
# TLEs = ST.read_tle_from_string(line1, line2)
