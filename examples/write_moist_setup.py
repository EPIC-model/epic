#!/usr/bin/env python
import numpy as np
from tools.nc_parcels import nc_parcels
import argparse

try:
    parser = argparse.ArgumentParser(
        description="Write parcels for moist setup."
    )

    parser.add_argument(
        "--ngrid",
        type=int,
        required=False,
        default=32,
        help="number of grid cells per dimension",
    )

    args = parser.parse_args()

    # should eventually read these from namelist?
    RH = 0.8
    z_c = 2500.0
    mu = 0.9
    z_d = 4000.0
    z_m = 5000.0
    r_plume = 800.0
    e_values = np.array([0.3, -0.4, 0.5])
    l_condense = 1000.0
    q0 = 0.015
    theta_l0 = 288.0
    ngrid = args.ngrid
    n_par_res = 2  # how many parcels per grid box in each dimension
    gravity = 9.81
    L_v = 2.501e6
    c_p = 1005.0

    parcels_per_dim = ngrid * n_par_res
    tuple_origin = (0, 0, 0)
    tuple_extent = (6280., 6280., 6280.)
    tuple_ncells = (np.int32(ngrid), np.int32(ngrid), np.int32(ngrid))
    origin=np.array(tuple_origin)
    extent=np.array(tuple_extent)
    ncells=np.array(tuple_ncells)

    dx_parcel = extent / parcels_per_dim

    ncp = nc_parcels()

    ncp.open('moist_parcels_restart.nc')

    ncp.add_box(tuple_origin, tuple_extent, tuple_ncells)

    # generate data
    if RH > 1.0:
        print("Error: Relative humidity fraction must be < 1.")
        quit()

    h_pl = q0 * np.exp(-z_c / l_condense)

    print("Humidity inside the plume is " + str(h_pl))

    if (mu > 1.0) or (mu <= RH):
        print("Error: mu must be between H and 1. The selected value is " + str(mu))
        quit()

    h_bg = mu * h_pl

    print("Background humidity is " + str(h_bg))

    z_b = l_condense * np.log(q0 * RH / h_bg)

    print("Base of mixed layer is " + str(z_b))

    dbdz = (
        (gravity * L_v / (c_p * theta_l0))
        * (h_pl - q0 * np.exp(-z_m / l_condense))
        / (z_m - z_d)
    )

    print("The buoyancy frequency in the stratified zone is " + str(np.sqrt(dbdz)))

    # Also obtain the plume liquid-water buoyancy (using also z_b):
    b_pl = dbdz * (z_d - z_b)

    print("The plume liquid water buoyancy b_pl = " + str(b_pl))
    print("corresponding to (theta_l-theta_l0)/theta_l0 = " + str(b_pl * gravity))

    if 2.0 * r_plume > z_b:
        print("Error: Plume radius is too big. At most it can be " + str(0.5 * z_b))
        quit()

    radsq = r_plume ** 2
    e_values = e_values / radsq

    if origin[2] > 0.0:
        print("The vertical origin must be zero.")
        quit()

    print("Box layout:")
    print("z_max           =" + str(extent[2]))
    print("z_m             =" + str(z_m))
    print("z_d             =" + str(z_d))
    print("z_c             =" + str(z_c))
    print("z_b             =" + str(z_b))
    print("zmin            =" + str(origin[2]))
    print("top of plume    =" + str(2.0 * r_plume))
    print("bottom of plume =" + str(0.0))

    centre = 0.5 * (2.0 * origin + extent)
    num_parcel = parcels_per_dim ** 3
    position = np.zeros((num_parcel, 3))
    buoyancy = np.zeros(num_parcel)
    humidity = np.zeros(num_parcel)
    volume = np.ones(num_parcel) * dx_parcel[0] * dx_parcel[1] * dx_parcel[2]
    b_diag = (0.75 * dx_parcel[0] * dx_parcel[1] * dx_parcel[2] / np.pi) ** (2.0 / 3.0)

    B = np.zeros((num_parcel, 5))
    B[:, 0] = b_diag
    B[:, 3] = b_diag
    vorticity = np.zeros((num_parcel, 3))

    iparcel = 0
    for ii in range(parcels_per_dim):
        for jj in range(parcels_per_dim):
            for kk in range(parcels_per_dim):
                pos = origin + dx_parcel * np.array([ii + 0.5, jj + 0.5, kk + 0.5])
                position[iparcel, :] = pos
                rpos1 = pos[0] - centre[0]
                rpos2 = pos[1] - centre[1]
                rpos3 = pos[2] - r_plume
                r2 = rpos1 ** 2 + rpos2 ** 2 + rpos3 ** 2
                if r2 <= radsq:
                    buoyancy[iparcel] = b_pl * (
                        1.0
                        + e_values[0] * rpos1 * rpos2
                        + e_values[1] * rpos1 * rpos3
                        + e_values[2] * rpos2 * rpos3
                    )
                    humidity[iparcel] = h_pl
                else:
                    if pos[2] < z_b:
                        # Mixed layer
                        buoyancy[iparcel] = 0.0
                        humidity[iparcel] = h_bg
                    else:
                        # Stratified layer
                        buoyancy[iparcel] = dbdz * (pos[2] - z_b)
                        humidity[iparcel] = q0 * RH * np.exp(-pos[2] / l_condense)
                iparcel = iparcel + 1

    # write parcel attributes
    ncp.add_dataset('x_position', position[:, 0], unit='m')
    ncp.add_dataset('y_position', position[:, 1], unit='m')
    ncp.add_dataset('z_position', position[:, 2], unit='m')

    ncp.add_dataset('buoyancy', buoyancy, unit='m/s^2')

    ncp.add_dataset('humidity', humidity, unit='1')

    ncp.add_dataset('volume', volume, unit='m^3')

    ncp.add_dataset('x_vorticity', vorticity[:, 0], unit='1/s')
    ncp.add_dataset('y_vorticity', vorticity[:, 1], unit='1/s')
    ncp.add_dataset('z_vorticity', vorticity[:, 2], unit='1/s')

    ncp.add_dataset('B11', B[:, 0], unit='m^2')
    ncp.add_dataset('B12', B[:, 1], unit='m^2')
    ncp.add_dataset('B13', B[:, 2], unit='m^2')
    ncp.add_dataset('B22', B[:, 3], unit='m^2')
    ncp.add_dataset('B23', B[:, 4], unit='m^2')

    ncp.add_physical_constant('q_0', q0)
    ncp.add_physical_constant('liquid_water_potential_temperature', theta_l0)
    ncp.add_physical_constant('gravity', gravity)
    ncp.add_physical_constant('latent_heat', L_v)
    ncp.add_physical_constant('specific_heat', c_p)

    ncp.close()
except Exception as ex:
    print(ex)
