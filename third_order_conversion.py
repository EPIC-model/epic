from numba import njit, prange, set_num_threads
import netCDF4 as nc
import numpy as np
from glob import glob
from math import erf
import argparse
import os
# Note: check division by zero errors.

projfac =8  # Subgrid zoom factor
r_limit_fac = 3.0
#r_limit_fac = 1.0 # Alternative: individual ellipsoid visualisation
len_condense = 1000.0
q_scale = 0.015
nthreads = 4
set_num_threads(nthreads)

parser = argparse.ArgumentParser()

parser.add_argument("input_file_name", type=str)
parser.add_argument("time_step", type=int)
args = parser.parse_args()
input_file_name = args.input_file_name
time_step = args.time_step

file_root, file_ext = os.path.splitext(input_file_name)

if not (file_ext == ".hdf5"):
    raise argparse.ArgumentTypeError('Argument filename must end on ".hdf5"')

ds_nc = nc.Dataset(input_file_name)

# set up time array
tstep_keys = []
for key in ds_nc.groups.keys():
    if "step#" in key:
        tstep_keys.append(key)

ds_ts = ds_nc.groups[tstep_keys[time_step]]

# set up spatial arrays
extent = ds_nc.groups["box"].extent
origin = ds_nc.groups["box"].origin
ncells = ds_nc.groups["box"].ncells

# initialisation
nx = ncells[0]
ny = ncells[1]
nz = ncells[2]  # Don't add 1 (this is non-periodic) here
# projection grid
nxproj = nx * projfac
nyproj = ny * projfac
nzproj = nz * projfac + 1  # But do add 1 here
dx = extent[0] / nx
dy = extent[1] / ny
dz = extent[2] / nz
dx_project = dx / projfac
dy_project = dy / projfac
dz_project = dz / projfac
dxi_project = projfac / dx  # Inverse of projection grid distance
dyi_project = projfac / dy
dzi_project = projfac / dz

# Projection grid coordinates
xp = np.linspace(origin[0], origin[0] + extent[0], nxproj + 1)[0:nxproj]
yp = np.linspace(origin[1], origin[1] + extent[1], nyproj + 1)[0:nyproj]
zp = np.linspace(origin[2], origin[2] + extent[2], nzproj)

ipi = 1.0 / np.pi
# Summation fields
xyz_field = np.zeros((nxproj, nyproj, nzproj), np.double)
xyz_volg = np.zeros((nxproj, nyproj, nzproj), np.double)  # Normalisation factor

# This function operates on small subdomains
# Does not yet work in parallel
@njit(parallel=True)
def add_data(
    h,
    x,
    y,
    z,
    bb1,
    bb2,
    bb3,
    bb4,
    bb5,
    vol,
    xyz_field_thread,
    xyz_volg_thread,
    xlowerfile,
    ylowerfile,
):
    # Construction below implements something that is a lot like OpenMP reduction in numba
    # Each thread adds to its own field, and only later we add threads together in add_fields
    for i_thread in prange(nthreads):
        bb = np.zeros((3, 3), np.double)
        bbinv = np.zeros((3, 3), np.double)
        xarr = np.zeros((3, 1), np.double)
        evals = np.zeros(3, np.double)
        for i in range(len(h)):
            if (
                not i % nthreads == i_thread
            ):  # Use modulus here, not sure if contiguous blocks could be better
                continue
            vol_now = vol[i]
            r_now = (vol_now * 0.75 * ipi) ** (1.0 / 3.0)
            # Maximum radius for parcel projection, and related constants
            x_now = x[i]
            y_now = y[i]
            z_now = z[i]
            scalar_now = h[i]
            b1 = bb1[i]
            b2 = bb2[i]
            b3 = bb3[i]
            b4 = bb4[i]
            b5 = bb5[i]
            abc = 0.75 * vol_now * ipi
            b6 = (
                abc * abc - b3 * (b2 * b5 - b3 * b4) + b1 * b5 * b5 - b2 * b3 * b5
            ) / (b1 * b4 - b2 * b2)
            bb[0, 0] = b1
            bb[0, 1] = b2
            bb[0, 2] = b3
            bb[1, 0] = b2
            bb[1, 1] = b4
            bb[1, 2] = b5
            bb[2, 0] = b3
            bb[2, 1] = b5
            bb[2, 2] = b6
            evals[:] = np.linalg.eigvalsh(bb)
            bbinv[:, :] = np.linalg.inv(bb)
            anisotropy_fact = np.sqrt(np.nanmax(evals)) / (
                (evals[0] * evals[1] * evals[2]) ** (1.0 / 6.0)
            )
            r_max_safety = r_limit_fac * anisotropy_fact
            r_max = r_max_safety * r_now
            # Determine which part of the grid the parcel contributes to
            xlower = int(dxi_project * ((x_now - r_max) - origin[0])) + 1
            xupper = int(dxi_project * ((x_now + r_max) - origin[0]))
            ylower = int(dyi_project * ((y_now - r_max) - origin[1])) + 1
            yupper = int(dyi_project * ((y_now + r_max) - origin[1]))
            zlower = max(int(dzi_project * ((z_now - r_max) - origin[2])), 0)
            zupper = min(int(dzi_project * ((z_now + r_max) - origin[2])), nzproj - 1)
            norm_fact = 1.0
            # Include upper values in loops (converted from Fortran)
            for ii in range(xlower, xupper + 1):
                for jj in range(ylower, yupper + 1):
                    for kk in range(zlower, zupper + 1):
                        xdist = ii * dx_project + origin[0] - x_now
                        ydist = jj * dy_project + origin[1] - y_now
                        zdist = kk * dz_project + origin[2] - z_now
                        dist2 = xdist * xdist + ydist * ydist + zdist * zdist
                        if dist2 < r_max * r_max:
                            xarr[0, 0] = xdist
                            xarr[1, 0] = ydist
                            xarr[2, 0] = zdist
                            r_act = np.sqrt(
                                np.dot(
                                    np.dot(xarr.transpose(), bbinv),
                                    xarr,
                                )[0, 0]
                            )
                            # r_act = np.sqrt(dist2)/r_now
                            if r_act < r_limit_fac:
                                # Correct y and x location with respect to file-based grids
                                # Do not correct for cyclic boundary conditions yet, this is done in add_fields
                                ii_corr = ii - xlowerfile
                                jj_corr = jj - ylowerfile
                                r_act3 = r_act * r_act * r_act
                                #r_fact = (r_act<1.0) # Alternative: individual ellipsoid visulatisation.
                                r_fact = np.exp(-r_act3)
                                xyz_field_thread[i_thread, ii_corr, jj_corr, kk] += (
                                    norm_fact * r_fact * scalar_now
                                )
                                xyz_volg_thread[i_thread, ii_corr, jj_corr, kk] += (
                                    norm_fact * r_fact
                                )


# Sum up thread-specific sums on local coordinates to the correct global coordinates
@njit()
def add_fields(
    xyz_field,
    xyz_volg,
    xyz_field_thread_file,
    xyz_volg_thread_file,
    xlowerfile,
    ylowerfile,
    nx_file,
    ny_file,
):
    for i_thread in range(nthreads):
        # Parallelise this routine in x, rather than over threads, so that conflicts are avoided
        for ii in range(nx_file):
            for jj in range(ny_file):
                for kk in range(nzproj):
                    itarg = (ii + xlowerfile) % nxproj  # Periodic BCs
                    jtarg = (jj + ylowerfile) % nyproj  # Periodic BCs
                    xyz_field[itarg, jtarg, kk] = (
                        xyz_field[itarg, jtarg, kk]
                        + xyz_field_thread_file[i_thread, ii, jj, kk]
                    )
                    xyz_volg[itarg, jtarg, kk] = (
                        xyz_volg[itarg, jtarg, kk]
                        + xyz_volg_thread_file[i_thread, ii, jj, kk]
                    )


h = ds_ts.variables["humidity"][:]
x = ds_ts.variables["position"][:, 0]
y = ds_ts.variables["position"][:, 1]
z = ds_ts.variables["position"][:, 2]
vol = ds_ts.variables["volume"][:]
bb1 = ds_ts.variables["B"][:, 0]
bb2 = ds_ts.variables["B"][:, 1]
bb3 = ds_ts.variables["B"][:, 2]
bb4 = ds_ts.variables["B"][:, 3]
bb5 = ds_ts.variables["B"][:, 4]

# Local domain min and max values
xminval = np.nanmin(x)
xmaxval = np.nanmax(x)
yminval = np.nanmin(y)
ymaxval = np.nanmax(y)
# Local domain projection coordinates boundaries
# Add one (round up rather than down), as values are zero at r_max distance
xlowerfile = int(dxi_project * ((xminval - 3.0 * dx) - origin[0])) + 1
xupperfile = int(dxi_project * ((xmaxval + 3.0 * dx) - origin[0]))
ylowerfile = int(dyi_project * ((yminval - 3.0 * dy) - origin[1])) + 1
yupperfile = int(dyi_project * ((ymaxval + 3.0 * dy) - origin[1]))
# Size of local projection domain, which includes both xupper and xlower etc.
nx_file = xupperfile - xlowerfile + 1
ny_file = yupperfile - ylowerfile + 1
# Initialise "thread-based" fields for this file
xyz_field_thread_file = np.zeros((nthreads, nx_file, ny_file, nzproj), np.double)
xyz_volg_thread_file = np.zeros((nthreads, nx_file, ny_file, nzproj), np.double)
add_data(
    h,
    x,
    y,
    z,
    bb1,
    bb2,
    bb3,
    bb4,
    bb5,
    vol,
    xyz_field_thread_file,
    xyz_volg_thread_file,
    xlowerfile,
    ylowerfile,
)
add_fields(
    xyz_field,
    xyz_volg,
    xyz_field_thread_file,
    xyz_volg_thread_file,
    xlowerfile,
    ylowerfile,
    nx_file,
    ny_file,
)

with open("tracker.txt", "a") as f:
    print("added fields", file=f)

#del xyz_field_thread_file, xyz_volg_thread_file

# Final calculation of liquid water field
# Normalize h field, and convert to LWC
# See MONC/MPIC intercomparison appendix B
# This is like strategy C
@njit()
def final_calc(field, volg, zp):
    for ii in range(np.shape(field)[0]):
        for jj in range(np.shape(field)[1]):
            for kk in range(np.shape(field)[2]):
                field_temp = field[ii, jj, kk] / volg[ii, jj, kk]
                volg_temp = max(
                    field_temp - q_scale * np.exp((origin[2] - zp[kk]) / len_condense),
                0.0,
                )
                field[ii, jj, kk] = field_temp
                volg[ii, jj, kk] = volg_temp


# store total water in field
# store liquid water in volg
final_calc(xyz_field, xyz_volg, zp)

with open("tracker.txt", "a") as f:
    print("div ready", file=f)

# Write to netcdf file
def write_field(field, field_name, out_file, mode):
    ncfile = nc.Dataset(out_file, mode, format="NETCDF4")
    if mode=="w":
        ncfile.createDimension("xp", nxproj)
        ncfile.createDimension("yp", nyproj)
        ncfile.createDimension("zp", nzproj)
        nc_xp = ncfile.createVariable("xp", "f8", ("xp",))
        nc_yp = ncfile.createVariable("yp", "f8", ("yp",))
        nc_zp = ncfile.createVariable("zp", "f8", ("zp",))
        nc_xp.units = "-"
        nc_yp.units = "-"
        nc_zp.units = "-"
        nc_xp[:] = xp
        nc_yp[:] = yp
        nc_zp[:] = zp
    nc_hl = ncfile.createVariable(field_name, np.dtype("f4").char, ("xp", "yp", "zp"))
    nc_hl.units = "-"
    nc_hl[:, :, :] = field[:, :, :]
    ncfile.close()

write_field(xyz_field, "hh", file_root + "_s3_" + str(time_step) + ".nc","w")
write_field(xyz_volg, "hl", file_root + "_s3_" + str(time_step) + ".nc","a")
