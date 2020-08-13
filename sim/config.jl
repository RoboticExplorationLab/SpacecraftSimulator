using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox
using YAML
const SD = SatelliteDynamics
const ST = SatelliteToolbox

function config(path_to_config_yaml)
      """This loads in the YAML for the sim run"""

      # load YAML
      data = YAML.load(open(joinpath(dirname(dirname(@__FILE__)),path_to_config_yaml)))

      # break it into the sections
      sc_properties = data["sc_properties"]
      gravity = data["gravity"]
      environment = data["environment"]
      init_oe = data["initial_orbital_elements"]
      propagator_parameters = data["propagator_parameters"]
      initial_attitude_conditions = data["initial_attitude_conditions"]
      inertia_properties = data["inertia_properties"]
      sensor_properties = data["sensor_properties"]

      # inertia matrix
      J = [sc_properties["inertia_rows"][1]';
           sc_properties["inertia_rows"][2]';
           sc_properties["inertia_rows"][3]']

      # store inverse of inertia matrix
      invJ = inv(J)

      # spacecraft faces
      faces = [sc_properties["faces"][1]';
               sc_properties["faces"][2]';
               sc_properties["faces"][3]';
               sc_properties["faces"][4]';
               sc_properties["faces"][5]';
               sc_properties["faces"][6]']

      # spacecraft properties
      sc_mass = sc_properties["mass"]
      sc_area = sc_properties["area"]
      sc_cd = sc_properties["cd"]

      # max magnetorquer dipole moments
      max_dipoles = sc_properties["max_dipoles"]

      # store all spacecraft properties in named tuple 'sc'
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
      t_vec_attitude = 0:dt_attitude:(t_vec_orbital[end])

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

      # initial orbital elements
      oe0 = [sma;ecc;inc;raan;argp;M]

      # convert to ECI r and v
      eci_rv_0 = sOSCtoCART(oe0, use_degrees=false)

      initial_conditions = (epc_orbital = epc_orbital,ᴺqᴮ0=ᴺqᴮ0,ω0=ω0,
                            eci_rv_0 = eci_rv_0,oe0 = oe0)

      # sample inertia
      if inertia_properties["sample_new"]
            J_sample = sample_inertia(J,
                                      inertia_properties["rotate_std_deg"],
                                      inertia_properties["scale_std"])
      else
            J_sample = copy(J)
      end

      # generate sensor offests
      gyro_properties = sensor_properties["gyro"]
      sun_sensor_properties = sensor_properties["sun_sensor"]
      magnetometer_properties = sensor_properties["magnetometer"]

      gyro = (rotate_std_deg   = gyro_properties["rotate_std_deg"],
              noise_std_degps  = gyro_properties["noise_std_degps"],
              initial_bias_deg = gyro_properties["initial_bias_deg"],
              bias_noise_std   = gyro_properties["bias_noise_std"])

      sun_sensor = (rotate_std_deg = sun_sensor_properties["rotate_std_deg"],
                    noise_std_deg  = sun_sensor_properties["noise_std_deg"])

      magnetometer=(rotate_std_deg = magnetometer_properties["rotate_std_deg"],
                    noise_std_deg  = magnetometer_properties["noise_std_deg"])

      offsets = (gyro = dcm_from_phi(deg2rad(gyro.rotate_std_deg)*randn(3)),
      gyro_bias = deg2rad(gyro.initial_bias_deg)*normalize(randn(3)),
      sun_sensor = dcm_from_phi(deg2rad(sun_sensor.rotate_std_deg)*randn(3)),
      magnetometer=dcm_from_phi(deg2rad(magnetometer.rotate_std_deg)*randn(3)))

      sensors = (gyro = gyro, sun_sensor = sun_sensor, magnetometer=magnetometer,
      offsets = offsets,J = J_sample, invJ = inv(J_sample))

      global params = (sc=sc, time_params=time_params,initial_conditions=initial_conditions,
                grav_deg = grav_deg, grav_order = grav_order,
                spherical_harmonic_gravity_bool = spherical_harmonic_gravity_bool,
                sensors = sensors)


      return initial_conditions, time_params
end



# read the TLE's
# line1 = "1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991"
# line2 = "2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102"
#
#
# TLEs = ST.read_tle_from_string(line1, line2)
