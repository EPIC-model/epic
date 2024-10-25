#!/usr/bin/env python
from tools.nc_fields import nc_fields
import numpy as np
import argparse

#The surface flux problem has few parameters, and it is not worth trying to make an exact comparison with Sam's results because he made some odd choices for some of the parameters.  The cleanest thing you can do is start with b = z and omega = 0, with f = 0 (rotation would be interesting, later).  Then apply a surface flux of 1- units do not matter and db/dz and the flux set the time and length scales to be O(1).  So, make the domain big enough, maybe 4L x 4L x L (in x, y and z) with L > 1, maybe 5.  Use a uniform flux but slightly disturb the initial parcel positions - this is enough.  I'll ask Sam to send his thesis to us.

try:
    parser = argparse.ArgumentParser(
        description="Create surface flux setup"
    )

    parser.add_argument(
        "--nx",
        type=int,
        required=False,
        default=256,
        help="number of cells in x",
    )

    parser.add_argument(
        "--ny",
        type=int,
        required=False,
        default=256,
        help="number of cells in y",
    )

    parser.add_argument(
        "--nz",
        type=int,
        required=False,
        default=64,
        help="number of cells in z",
    )

    parser.add_argument(
        "--L",
        type=float,
        required=False,
        default=3200.0,
        help="domain extent",
    )

    args = parser.parse_args()

    lapse=0.003
    thetafluxmag=0.1

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz

    L = args.L

    ncf = nc_fields()

    ncf.open('boundary_layer_setup' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.nc')

    # domain origin
    origin = (-2.0 * L, -2.0 * L, 0.0)

    # domain extent
    extent = (4.0 * L, 4.0 * L, L)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    buoy = np.zeros((nz+1, ny, nx))

    # ranges from 0 to nz
    for k in range(nz+1):
        zrel = k * dz
        buoy[k, :, :] = 300.+lapse*zrel

    # add perturbation
    rng = np.random.default_rng(seed=42)
    noise = rng.uniform(low=0.0, high=0.01*(buoy.max()-300.), size=(nz+1, ny, nx))
    buoy += noise

    # write all provided fields
    ncf.add_field('theta', buoy, unit='K', long_name='potential temperature')

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.close()

    #
    # Setup surface fluxes:
    #

    ncf_flux = nc_fields()

    ncf_flux.set_dim_names(['x', 'y'])

    ncf_flux.open('surface_flux_' + str(nx) + 'x' + str(ny) + '.nc')

    # domain origin
    origin = (-2.0 * L, -2.0 * L)

    # domain extent
    extent = (4.0 * L, 4.0 * L)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny

    thetaflux = np.ones((ny, nx))*thetafluxmag
    qvflux = np.zeros((ny, nx))

    ncf_flux.add_field('thetaflux', thetaflux, unit='K m/s')
    ncf_flux.add_field('qvflux', qvflux, unit='kg/kg m/s')

    ncf_flux.add_box(origin, extent, [nx, ny])

    ncf_flux.close()

except Exception as err:
    print(err)

