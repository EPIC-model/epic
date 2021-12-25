# transpose the fields too!
import xarray as xr
import argparse
import os.path
import netCDF4 as nc
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("input_file_name", type=str)
args = parser.parse_args()
input_file_name = args.input_file_name

file_root, file_ext = os.path.splitext(input_file_name)

if not (file_ext == ".hdf5"):
    raise argparse.ArgumentTypeError('Argument filename must end on ".hdf5"')


ds_nc=nc.Dataset(input_file_name)

# set up time array
tstep_keys=[]
for key in ds_nc.groups.keys():
    if 'step#' in key:
        tstep_keys.append(key)

tlen=len(tstep_keys)
time_array=np.zeros(tlen)
for i_time in range(tlen):
   time_array[i_time]=ds_nc.groups[tstep_keys[i_time]].t

# set up spatial arrays
extent=ds_nc.groups['box'].extent
origin=ds_nc.groups['box'].origin
ncells=ds_nc.groups['box'].ncells

x_array=np.linspace(origin[0],origin[0]+extent[0],ncells[0]+1)[:-1]
y_array=np.linspace(origin[1],origin[1]+extent[1],ncells[1]+1)[:-1]
z_array=np.linspace(origin[2],origin[2]+extent[2],ncells[2]+1)
elem_array=range(3)

step0=ds_nc.groups[tstep_keys[0]]

ds_xr = xr.Dataset()

for variable in step0.variables:
    n_dims=len(np.shape(step0[variable][:]))
    if(n_dims==3):
        ds_xr[variable] = xr.DataArray(coords=[time_array, z_array, y_array, x_array], dims=["time", "z","y","x"])
    elif(n_dims==4):
        ds_xr[variable] = xr.DataArray(coords=[time_array, z_array, y_array, x_array, elem_array], dims=["time", "z","y","x", "index"])

ds_xr['diffB'] = xr.DataArray(coords=[time_array, z_array, y_array, x_array], dims=["time", "z","y","x"])

for i_time in range(tlen):
   this_step=ds_nc.groups[tstep_keys[i_time]]
   for variable in this_step.variables:
        ds_xr[variable][i_time]=np.transpose(this_step[variable][:])
   ds_xr['diffB'][i_time]=np.transpose(this_step['total buoyancy'][:])-np.transpose(this_step['dry buoyancy'][:])

ds_xr.to_netcdf(file_root + ".nc")
