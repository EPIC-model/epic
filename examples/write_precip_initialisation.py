#!/usr/bin/env python
import numpy as np
import xarray as xr

r_plume = 800.0
centre = [0.0, 3000.0]
dx_parcel = 5.0 # Distance between parcels

n_steps = int(np.floor(r_plume / dx_parcel))
parcel_shifts = np.arange(-n_steps, n_steps + 1) * dx_parcel
n_shifts=len(parcel_shifts)

i_parcel=0
# oversized_initialisation
x_array=np.zeros((1,n_shifts*n_shifts))
z_array=np.zeros((1,n_shifts*n_shifts))

qr_parcels=0.002
Nr_parcels=10000
parcel_volume=dx_parcel*dx_parcel

for ii in range(n_shifts):
    parcel_shift_x=parcel_shifts[ii]
    for jj in range(n_shifts):
        parcel_shift_z=parcel_shifts[jj]
        if(parcel_shift_x*parcel_shift_x+parcel_shift_z*parcel_shift_z<r_plume*r_plume):
             x_array[0,i_parcel] = centre[0] + parcel_shift_x
             z_array[0,i_parcel] = centre[1] + parcel_shift_z
             i_parcel=i_parcel+1

x_array=x_array[:,:i_parcel]
z_array=z_array[:,:i_parcel]

len_parcels=np.shape(x_array)[1]

volume_array=np.ones((1,len_parcels))*parcel_volume
qr_array=np.ones((1,len_parcels))*qr_parcels
Nr_array=np.ones((1,len_parcels))*Nr_parcels

n_parcels = np.arange(1, len_parcels+1,dtype=np.int32)

time = np.array([0.0])

# Coordinates
coords = {
    "time": ("time", time, {
        "units": "seconds since 1970-01-01 00:00:00",
        "calendar": "proleptic_gregorian"
    }),
    "n_parcels": ("n_parcels", n_parcels)
}

# Create the actual dataset
ds = xr.Dataset(
    {
        "x_position": xr.DataArray(x_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m", "long_name": "x position component"}),
        "z_position": xr.DataArray(z_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m", "long_name": "z position component"}),
        "volume": xr.DataArray(volume_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "m^2", "long_name": "parcel volume"}),
        "qr": xr.DataArray(qr_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "kg/kg", "long_name": "rain mixing ratio"}),
        "Nr": xr.DataArray(Nr_array, dims=["time", "n_parcels"], coords=coords,
                              attrs={"units": "1/kg", "long_name": "rain number concentration"})
    }
)

# Save with unlimited time dimension
ds.to_netcdf("rain_initiation_dataset.nc", unlimited_dims=["time"])
