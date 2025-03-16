#!/usr/bin/env python
#
# Rayleigh-Taylor flow with buoyancy
#
#       b(x, y, z) = = -sin(z) + epsilon * cos^2(z) * h(x,y)
#
# and horizontal perturbation
#
#       h(x, y) = cos(4x) * cos(2y + pi/6) + sin(2x + pi/6) * sin(4y)
#
# in the domain [-pi/2, pi/2]^3.
from tools.nc_fields import nc_fields
import numpy as np
import argparse
import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
from tools.mpl_style import *

try:
    parser = argparse.ArgumentParser(
        description="Create Rayleigh-Taylor flow"
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
        default=64,
        help="number of cells in z",
    )

    parser.add_argument(
        "--epsilon",
        type=float,
        required=False,
        default=0.1,
        help="perturbation amplitude",
    )

    parser.add_argument(
        "--plot",
        type=bool,
        required=False,
        default=False,
        help='plot the horizontal perturbation'
    )

    parser.add_argument(
        "--ape-calculation",
        type=str,
        default='none',
        choices=['none', 'sorting', 'ape density'],
        help="Option to calculate the available potential energy; " + \
             "ensure src/utils/ape_density.f90 is implemented when enabled"
    )

    args = parser.parse_args()

    # number of cells
    nx = args.nx
    ny = args.ny
    nz = args.nz


    # perturbation amplitude
    eps = args.epsilon

    ncf = nc_fields()

    ncf.open('rt_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.nc')

    # domain origin
    origin = (-0.5 * np.pi, -0.5 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (np.pi, np.pi, np.pi)


    # mesh spacings
    dx = extent[0] / nx
    dy = extent[1] / ny
    dz = extent[2] / nz

    buoy = np.zeros((nz+1, ny, nx))

    h = np.zeros((ny, nx))

    # ranges from 0 to nx-1
    for i in range(nx):
        for j in range(ny):
            x = origin[0] + i * dx
            y = origin[1] + j * dy

            h[j, i] = np.cos(4.0*x) * np.cos(2.0*y + np.pi/6.0) + np.sin(2.0*x + np.pi/6.0) * np.sin(4.0*y)

    # ranges from 0 to nz
    for k in range(nz+1):
        z = origin[2] + k * dz
        buoy[k, :, :] = -np.sin(z) + eps * np.cos(z) ** 2 * h[:, :]

    # write all provided fields
    ncf.add_field('buoyancy', buoy, unit='m/s^2', long_name='buoyancy')

    ncf.add_box(origin, extent, [nx, ny, nz])

    ncf.add_axis('t', [0.0])

    ncf.add_physical_quantity('l_planetary_vorticity', 'true')
    ncf.add_physical_quantity('planetary_angular_velocity', 0.5)
    ncf.add_physical_quantity('latitude_degrees', 90)
    ncf.add_physical_quantity('ape_calculation', args.ape_calculation)

    ncf.close()

    if args.plot:
        mpl.rcParams['font.size'] = 16
        plt.figure(figsize=(7, 4), dpi=200)

        im = plt.imshow(h,
                        origin='lower',
                        interpolation='bilinear',
                        cmap=cc.cm['rainbow4'],
                        extent=[origin[0], origin[0]+extent[0],
                                origin[1], origin[1]+extent[1]])

        cbar = plt.colorbar(im, pad=0.03, ticks=[-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
        cbar.set_label(r'$h(x, y)$')

        ticks = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
        tlabs = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

        plt.xticks(ticks, tlabs)
        plt.yticks(ticks, tlabs)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')

        plt.savefig('horizontal_perturbation.pdf', bbox_inches='tight')
        plt.close()

except Exception as err:
    print(err)
