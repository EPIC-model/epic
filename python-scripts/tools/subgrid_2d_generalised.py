from numba import njit, prange, set_num_threads
import h5py
from tools.h5_reader import H5Reader
import numpy as np
from glob import glob
from math import erf
import argparse
import os

# Note: check division by zero errors.

projfac = 12  # Subgrid zoom factor
r_limit_fac = 4.0  # How far out we project the parcels, expressed in ellipsoid radii.
nthreads = 1
set_num_threads(nthreads)

parser = argparse.ArgumentParser()

parser.add_argument("input_file_name", type=str)
parser.add_argument("a_fac", type=float)
parser.add_argument("b_fac", type=float)
parser.add_argument("time_step", type=int)

args = parser.parse_args()
input_file_name = args.input_file_name
time_step = args.time_step
a_fac = args.a_fac
b_fac = args.b_fac

file_root, file_ext = os.path.splitext(input_file_name)

if not (file_ext == ".hdf5"):
    raise argparse.ArgumentTypeError('Argument filename must end on ".hdf5"')

h5reader = H5Reader()
h5reader.open(input_file_name)


# set up spatial arrays
extent = h5reader.get_box_extent()
origin = h5reader.get_box_origin()
ncells = h5reader.get_box_ncells()

# initialisation
nx = ncells[0]
nz = ncells[1]  # Don't add 1 (this is non-periodic) here
# projection grid
nxproj = nx * projfac
nzproj = nz * projfac + 1  # But do add 1 here
dx = extent[0] / nx
dz = extent[1] / nz
dx_project = dx / projfac
dz_project = dz / projfac
dxi_project = projfac / dx  # Inverse of projection grid distance
dzi_project = projfac / dz

# Projection grid coordinates
xp = np.linspace(origin[0], origin[0] + extent[0], nxproj + 1)[0:nxproj]
zp = np.linspace(origin[1], origin[1] + extent[1], nzproj)

ipi = 1.0 / np.pi
# Summation fields
xz_field = np.zeros((nxproj, nzproj), np.double)
xz_volg = np.zeros((nxproj, nzproj), np.double)  # Normalisation factor

# This function operates on small subdomains
# Does not yet work in parallel
@njit(parallel=True)
def add_data(
    scalar, x, z, bb1, bb2, vol, xz_field_thread, xz_volg_thread, xlowerfile,
):
    # Construction below implements something that is a lot like OpenMP reduction in numba
    # Each thread adds to its own field, and only later we add threads together in add_fields
    for i_thread in prange(nthreads):
        bb = np.zeros((2, 2), np.double)
        bbinv = np.zeros((2, 2), np.double)
        xarr = np.zeros((2, 1), np.double)
        evals = np.zeros(2, np.double)
        for i in range(len(scalar)):
            if (
                not i % nthreads == i_thread
            ):  # Use modulus here, not sure if contiguous blocks could be better
                continue
            vol_now = vol[i]
            r_now = (vol_now * ipi) ** (1.0 / 2.0)
            # Maximum radius for parcel projection, and related constants
            x_now = x[i]
            z_now = z[i]
            scalar_now = scalar[i]
            b1 = bb1[i]
            b2 = bb2[i]
            ab = vol_now * ipi
            b3 = (ab ** 2 + b2 ** 2) / b1
            bb[0, 0] = b1
            bb[0, 1] = b2
            bb[1, 0] = b2
            bb[1, 1] = b3
            evals[:] = np.linalg.eigvalsh(bb)
            bbinv[:, :] = np.linalg.inv(bb)
            anisotropy_fact = np.sqrt(np.nanmax(evals)) / (
                (evals[0] * evals[1]) ** (1.0 / 4.0)
            )
            r_max_safety = r_limit_fac * anisotropy_fact
            r_max = r_max_safety * r_now
            # Determine which part of the grid the parcel contributes to
            xlower = int(dxi_project * ((x_now - r_max) - origin[0])) + 1
            xupper = int(dxi_project * ((x_now + r_max) - origin[0]))
            zlower = max(int(dzi_project * ((z_now - r_max) - origin[1])), 0)
            zupper = min(int(dzi_project * ((z_now + r_max) - origin[1])), nzproj - 1)
            norm_fact = 1.0
            # Include upper values in loops (converted from Fortran)
            for ii in range(xlower, xupper + 1):
                for kk in range(zlower, zupper + 1):
                    xdist = ii * dx_project + origin[0] - x_now
                    zdist = kk * dz_project + origin[1] - z_now
                    dist2 = xdist * xdist + zdist * zdist
                    if dist2 < r_max * r_max:
                        xarr[0, 0] = xdist
                        xarr[1, 0] = zdist
                        r_act = np.sqrt(
                            np.dot(np.dot(xarr.transpose(), bbinv), xarr,)[0, 0]
                        )
                        # r_act = np.sqrt(dist2)/r_now
                        if r_act < r_limit_fac:
                            # Correct y and x location with respect to file-based grids
                            # Do not correct for cyclic boundary conditions yet, this is done in add_fields
                            ii_corr = ii - xlowerfile
                            # r_fact = (
                            #    r_act < 1.0
                            # )  # Alternative: individual ellipsoid visulatisation.
                            r_exponent = (r_act / a_fac) ** b_fac
                            r_fact = np.exp(-r_exponent)
                            xz_field_thread[i_thread, ii_corr, kk] += (
                                norm_fact * r_fact * scalar_now
                            )
                            xz_volg_thread[i_thread, ii_corr, kk] += norm_fact * r_fact


# Sum up thread-specific sums on local coordinates to the correct global coordinates
@njit()
def add_fields(
    xz_field, xz_volg, xz_field_thread_file, xz_volg_thread_file, xlowerfile, nx_file,
):
    for i_thread in range(nthreads):
        # Parallelise this routine in x, rather than over threads, so that conflicts are avoided
        for ii in range(nx_file):
            for kk in range(nzproj):
                itarg = (ii + xlowerfile) % nxproj  # Periodic BCs
                xz_field[itarg, kk] = (
                    xz_field[itarg, kk] + xz_field_thread_file[i_thread, ii, kk]
                )
                xz_volg[itarg, kk] = (
                    xz_volg[itarg, kk] + xz_volg_thread_file[i_thread, ii, kk]
                )


# Final calculation of liquid water field
# Normalize scalar field, and convert to LWC
# See MONC/MPIC intercomparison appendix B
# This is like strategy C
@njit()
def final_calc(field, volg, zp):
    for ii in range(np.shape(field)[0]):
        for kk in range(np.shape(field)[1]):
            if volg[ii, kk] > 0.0:
                field_temp = field[ii, kk] / volg[ii, kk]
            else:
                field_temp = np.nan
            field[ii, kk] = field_temp


# Write to netcdf file
def write_field(h5handle, field, field_name):
    h5handle[field_name] = np.float64(field.copy())

x = h5reader.get_dataset(time_step, "position")[:, 0]
z = h5reader.get_dataset(time_step, "position")[:, 1]
vol = h5reader.get_dataset(time_step, "volume")[:]
bb1 = h5reader.get_dataset(time_step, "B")[:, 0]
bb2 = h5reader.get_dataset(time_step, "B")[:, 1]

# Local domain min and max values
xminval = np.nanmin(x)
xmaxval = np.nanmax(x)
# Local domain projection coordinates boundaries
# Add one (round up rather than down), as values are zero at r_max distance
xlowerfile = int(dxi_project * ((xminval - 3.0 * dx) - origin[0])) + 1
xupperfile = int(dxi_project * ((xmaxval + 3.0 * dx) - origin[0]))
# Size of local projection domain, which includes both xupper and xlower etc.
nx_file = xupperfile - xlowerfile + 1



outfile = file_root + "_" + str(a_fac) + "_" + str(b_fac) + "_" + str(time_step) + "_fields.hdf5"

h5file = h5py.File(outfile, 'w')
dt = h5py.string_dtype('ascii', 6)
h5file.attrs.create("output_type", r"fields", dtype=dt, shape=1)
h5file.attrs['nsteps'] = [1]

box = h5file.create_group('box')
box.attrs['extent'] = extent
box.attrs['origin'] = origin
box.attrs['ncells'] = (np.int32(nxproj), np.int32(nzproj))

group = h5file.create_group('step#' + '0'.zfill(10))
group.attrs['t'] = 0.0


for target_variable in ["buoyancy"]: #, "vorticity"]:
    scalar = h5reader.get_dataset(time_step, target_variable).copy()

    # Initialise "thread-based" fields for this file
    xz_field_thread_file = np.zeros((nthreads, nx_file, nzproj), np.double)
    xz_volg_thread_file = np.zeros((nthreads, nx_file, nzproj), np.double)

    add_data(
        scalar,
        x,
        z,
        bb1,
        bb2,
        vol,
        xz_field_thread_file,
        xz_volg_thread_file,
        xlowerfile,
    )
    add_fields(
        xz_field,
        xz_volg,
        xz_field_thread_file,
        xz_volg_thread_file,
        xlowerfile,
        nx_file,
    )

    del xz_field_thread_file, xz_volg_thread_file

    final_calc(xz_field, xz_volg, zp)

    write_field(
        group,
        xz_field,
        target_variable
    )
#write_field(
    #group,
    #xz_volg,
    #"volume"
##)

h5file.close()

# python subgrid_2d_generalised.py straka_parcels.hdf5 1.0 2.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.2 2.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.5 2.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 2.0 2.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.0 2.5 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.2 2.5 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.5 2.5 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 2.0 2.5 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.0 3.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.2 3.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 1.5 3.0 45
# python subgrid_2d_generalised.py straka_parcels.hdf5 2.0 3.0 45
