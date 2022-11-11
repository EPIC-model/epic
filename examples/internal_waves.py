#!/usr/bin/env python
#

from tools.nc_fields import nc_fields
import numpy as np
import argparse

try:
    parser = argparse.ArgumentParser(
        description="Create flow"
    )

    parser.add_argument(
        "--nx",
        type=int,
        required=False,
        default=64,
        help="number of cells in x",
    )

    parser.add_argument(
        "--ny",
        type=int,
        required=False,
        default=64,
        help="number of cells in y",
    )

    parser.add_argument(
        "--nz",
        type=int,
        required=False,
        default=16,
        help="number of cells in z",
    )

    parser.add_argument(
        "--k",
        type=float,
        required=False,
        default=0.5,
        help="k wavenumber",
    )

    parser.add_argument(
        "--l",
        type=float,
        required=False,
        default=0.5,
        help="l wavenumber",
    )

    parser.add_argument(
        "--m",
        type=float,
        required=False,
        default=1,
        help="m wavenumber",
    )

    parser.add_argument(
        "--N",
        type=float,
        required=False,
        default=2.0,
        help="buoyancy frequency",
    )

    parser.add_argument(
        "--f",
        type=float,
        required=False,
        default=1.0,
        help="Coriolis frequency",
    )

    parser.add_argument(
        "--what",
        type=float,
        required=False,
        default=0.001,
        help="w amplitude",
    )

    args = parser.parse_args()

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz

    kk = args.k
    ll = args.l
    mm = args.m

    N = args.N
    f = args.f
    what = args.what

    N2 = N ** 2
    f2 = f ** 2
    sigma2 = (N2 * (kk ** 2 + ll ** 2) + f2 * mm ** 2) / (kk ** 2 + ll ** 2 + mm ** 2)

    ncf = nc_fields()

    ncf.open('internal_wave_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.nc')

    # domain origin
    origin = (-2.0 * np.pi, -2.0 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (4.0 * np.pi, 4.0 * np.pi, np.pi)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    xi = np.zeros((nz+1, ny, nx))
    eta = np.zeros((nz+1, ny, nx))
    zeta = np.zeros((nz+1, ny, nx))
    buoy = np.zeros((nz+1, ny, nx))


    N2f2 = N2 - f2
    s2f2i = 1.0 / (sigma2 - f2)
    sigma = np.sqrt(sigma2)
    sigi = 1.0 / sigma
    N2s2si = (N2 - sigma2) * sigi


    # ranges from 0 to nx-1
    for i in range(nx):
        for j in range(ny):
            for k in range(nz+1):
                x = origin[0] + i * dx
                y = origin[1] + j * dy
                z = origin[2] + k * dz

                phi = kk * x + ll * y

                sinphi = np.sin(phi)
                cosphi = np.cos(phi)
                cosmz = np.cos(mm * z)

                xi[k, j, i] = s2f2i * what * cosmz * (f * kk * N2s2si * cosphi - ll * N2f2 * sinphi)

                eta[k, j, i] = s2f2i * what * cosmz * (f * ll * N2s2si * cosphi + kk * N2f2 * sinphi)

                zeta[k, j, i] = f * mm * what * sigi * np.sin(mm * z) * sinphi

                buoy[k, j, i] = N2 * z + N2 * what * sigi * cosmz * sinphi


    # write all provided fields
    ncf.add_field('buoyancy', buoy, unit='m/s^2', long_name='buoyancy')
    ncf.add_field('x_vorticity', xi, unit='1/s', long_name='x-vorticity')
    ncf.add_field('y_vorticity', eta, unit='1/s', long_name='y-vorticity')
    ncf.add_field('z_vorticity', zeta, unit='1/s', long_name='z-vorticity')

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.add_physical_quantity('planetary_vorticity', 'true')
    ncf.add_physical_quantity('planetary_angular_velocity', 0.5 * f)
    ncf.add_physical_quantity('latitude_degrees', 90)

    ncf.close()

except Exception as err:
    print(err)
