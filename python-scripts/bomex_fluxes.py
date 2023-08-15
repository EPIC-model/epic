#!/usr/bin/env python
from tools.nc_fields import nc_fields
import numpy as np
import argparse

try:
    parser = argparse.ArgumentParser(description="Create surface flux setup")

    parser.add_argument(
        "--nx",
        type=int,
        required=False,
        default=192,
        help="number of cells in x",
    )

    parser.add_argument(
        "--ny",
        type=int,
        required=False,
        default=192,
        help="number of cells in y",
    )

    parser.add_argument(
        "--nz",
        type=int,
        required=False,
        default=96,
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

    thetafluxmag = 8.0e-3
    qvfluxmag = 5.2e-5

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz

    L = args.L

    ncf = nc_fields()

    ncf.open("bomex_setup" + str(nx) + "x" + str(ny) + "x" + str(nz) + ".nc")

    # domain origin
    origin = (-1.0 * L, -1.0 * L, 0.0)

    # domain extent
    extent = (2.0 * L, 2.0 * L, L)

    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    theta = np.zeros((nz + 1, ny, nx))
    qv = np.zeros((nz + 1, ny, nx))
    yvort = np.zeros((nz + 1, ny, nx))
    xvort = np.zeros((nz + 1, ny, nx))
    novort = np.zeros((nz + 1, ny, nx))

    # ranges from 0 to nz
    for k in range(nz + 1):
        zz = k * dz
        if zz < 520.0:
            theta[k, :, :] = 298.7
        elif zz < 1480:
            theta[k, :, :] = 298.7 + (zz - 520.0) * (302.4 - 298.7) / (1480.0 - 520.0)
        elif zz < 2000:
            theta[k, :, :] = 302.4 + (zz - 1480.0) * (308.2 - 302.4) / (2000.0 - 1480.0)
        else:
            theta[k, :, :] = 308.2 + (zz - 2000.0) * (311.85 - 308.2) / (
                3000.0 - 2000.0
            )

        # specific humidity
        if zz < 520.0:
            qv[k, :, :] = 1e-3 * (17.0 + zz * (16.3 - 17.0) / 520.0)
        elif zz < 1480:
            qv[k, :, :] = 1.0e-3 * (
                16.3 + (zz - 520.0) * (10.7 - 16.3) / (1480.0 - 520.0)
            )
        elif zz < 2000:
            qv[k, :, :] = 1.0e-3 * (
                10.7 + (zz - 1480.0) * (4.2 - 10.7) / (2000.0 - 1480.0)
            )
        else:
            qv[k, :, :] = 1.0e-3 * (
                4.2 + (zz - 2000.0) * (3.0 - 4.2) / (3000.0 - 2000.0)
            )

        if zz < 300.0:
            xvort[k, :, :] = 0.0043 * np.cos(np.pi * (zz + 600) / 1200.0) + 0.0023
            yvort[k, :, :] = 0.0062 * np.cos((np.pi * (zz - 300) / 600.0) ** 2) - 0.0067
        elif zz < 600.0:
            xvort[k, :, :] = 0.0043 * np.cos(np.pi * (zz + 600) / 1200.0) + 0.0023
            yvort[k, :, :] = 0.0027 * np.cos(np.pi * (zz - 300) / 600.0) - 0.0032
        elif zz < 1400.0:
            xvort[k, :, :] = -0.001 - 0.001 * np.cos(np.pi * (zz - 600.0) / 800.0)
            yvort[k, :, :] = 0.005 * np.cos(np.pi * (zz - 1400.0) / 1600.0) - 0.0032
        else:
            xvort[k, :, :] = 0.0
            yvort[k, :, :] = 0.0018

    # add perturbation
    rng = np.random.default_rng(seed=42)
    noise = rng.uniform(low=-0.05, high=0.05, size=(nz + 1, ny, nx))
    theta += noise

    # write all provided fields
    ncf.add_field("x_vorticity", xvort, unit="1/s")
    ncf.add_field("y_vorticity", yvort, unit="1/s")
    ncf.add_field("z_vorticity", novort, unit="1/s")
    ncf.add_field("theta", theta, unit="K", long_name="potential temperature")
    ncf.add_field("qv", qv, unit="kg/kg", long_name="water vapour spec. hum.")

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.close()

    #
    # Setup surface fluxes:
    #

    ncf_flux = nc_fields()

    ncf_flux.set_dim_names(["x", "y"])

    ncf_flux.open("bomex_flux_" + str(nx) + "x" + str(ny) + ".nc")

    # domain origin
    origin = (-1.0 * L, -1.0 * L)

    # domain extent
    extent = (2.0 * L, 2.0 * L)

    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny

    thetaflux = np.ones((ny, nx)) * thetafluxmag
    qvflux = np.ones((ny, nx)) * qvfluxmag

    ncf_flux.add_field("thetaflux", thetaflux, unit="K m/s")
    ncf_flux.add_field("qvflux", qvflux, unit="kg/kg m/s")

    ncf_flux.add_box(origin, extent, [nx, ny])

    ncf_flux.close()

except Exception as err:
    print(err)
