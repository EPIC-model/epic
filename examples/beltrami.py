#!/usr/bin/env python
#
# 3D Beltrami flow:
#   velocity:
#     u(x, y, z) = 1/4 * [sin(z) - 3 * cos(z)] * sin(2x + 2y)
#     v(x, y, z) = 1/4 * [sin(z) + 3 * cos(z)] * sin(2x + 2y)
#     w(x, y, z) = cos(z) * cos(2x + 2y)
#
#   vorticitiy:
#     \xi(x, y, z)   = 3 * u(x, y, z)
#     \eta(x, y, z)  = 3 * v(x, y, z)
#     \zeta(x, y, z) = 3 * w(x, y, z)
#
from tools.nc_fields import nc_fields
import numpy as np
import argparse

try:
    parser = argparse.ArgumentParser(
        description="Create Beltrami flow"
    )

    parser.add_argument(
        "--nx",
        type=int,
        required=False,
        default=16,
        help="number of cells in x",
    )

    parser.add_argument(
        "--ny",
        type=int,
        required=False,
        default=16,
        help="number of cells in y",
    )

    parser.add_argument(
        "--nz",
        type=int,
        required=False,
        default=16,
        help="number of cells in z",
    )

    args = parser.parse_args()

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz

    ncf = nc_fields()

    ncf.open('beltrami_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.nc')

    # domain origin
    origin = (-0.5 * np.pi, -0.5 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (np.pi, np.pi, np.pi)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    xi   = np.zeros((nz+1, ny, nx))
    eta  = np.zeros((nz+1, ny, nx))
    zeta = np.zeros((nz+1, ny, nx))

    # ranges from 0 to nx-1
    for i in range(nx):
        for j in range(ny):
            # ranges from 0 to nz
            for k in range(nz+1):
                x = origin[0] + i * dx
                y = origin[1] + j * dy
                z = origin[2] + k * dz

                # velocity:
                u = 0.25 * (np.sin(z) - 3.0 * np.cos(z)) * np.sin(2.0 * x + 2.0 * y)
                v = 0.25 * (np.sin(z) + 3.0 * np.cos(z)) * np.sin(2.0 * x + 2.0 * y)
                w = np.cos(z) * np.cos(2.0 * x + 2.0 * y)

                # vorticity:
                xi[k, j, i]   = 3.0 * u
                eta[k, j, i]  = 3.0 * v
                zeta[k, j, i] = 3.0 * w

    # write all provided fields
    ncf.add_field('x_vorticity', xi, unit='1/s')
    ncf.add_field('y_vorticity', eta, unit='1/s')
    ncf.add_field('z_vorticity', zeta, unit='1/s')

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.close()

except Exception as err:
    print(err)
