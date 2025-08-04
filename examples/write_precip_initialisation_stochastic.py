#!/usr/bin/env python
import numpy as np
import xarray as xr
from scipy.special import gamma
from scipy.stats import gamma as gamma_dist

# ---------- Functions ----------
# Initially written with GPT assistancei by SB
def lambda_from_nr_qr(nr, qr, rho_water, rho_air, mu):
    """Compute slope parameter lambda for CASIM gamma PSD."""
    return ((np.pi * nr * rho_water * (mu + 3) * (mu + 2) * (mu + 1)) /
            (6.0 * qr * rho_air )) ** (1.0/3.0)

def sample_mass_weighted(nr, qr, rho_water, rho_air, mu, size=1):
    """Draw diameters weighted by mass from CASIM gamma PSD."""
    lam = lambda_from_nr_qr(nr, qr, rho_water, rho_air, mu)
    shape_mass = mu + 4.0
    scale = 1.0 / lam
    return gamma_dist.rvs(a=shape_mass, scale=scale, size=size)

def equivalent_number_concentration(D, qr, rho_water, rho_air):
    """Equivalent N_r if all mass q were in droplets of diameter D."""
    mass_per_drop = (np.pi / 6.0) * rho_water * D**3
    return rho_air * qr / mass_per_drop

r_plume = 800.0
centre = [0.0, 3000.0]
dx_parcel = 10.0 # Distance between parcels

n_steps = int(np.floor(r_plume / dx_parcel))
parcel_shifts = np.arange(-n_steps, n_steps + 1) * dx_parcel
n_shifts=len(parcel_shifts)

qr_parcels=0.002
Nr_parcels=10000
mu=2.5
rho_water = 1000.0
rho_air = 1.2256 
multiplicity=10
parcel_volume=dx_parcel*dx_parcel

i_parcel=0
# oversized_initialisation
x_array=np.zeros((1,n_shifts*n_shifts*multiplicity))
z_array=np.zeros((1,n_shifts*n_shifts*multiplicity))

for ii in range(n_shifts):
    parcel_shift_x=parcel_shifts[ii]
    for jj in range(n_shifts):
        parcel_shift_z=parcel_shifts[jj]
        if(parcel_shift_x*parcel_shift_x+parcel_shift_z*parcel_shift_z<r_plume*r_plume):
             for mm in range(multiplicity):
                 x_array[0,i_parcel] = centre[0] + parcel_shift_x
                 z_array[0,i_parcel] = centre[1] + parcel_shift_z
                 i_parcel=i_parcel+1

x_array=x_array[:,:i_parcel]
z_array=z_array[:,:i_parcel]

len_parcels=np.shape(x_array)[1]

volume_array=np.ones((1,len_parcels))*parcel_volume
qr_array=np.ones((1,len_parcels))*qr_parcels/multiplicity
Nr_array=np.ones((1,len_parcels))

D_samples = sample_mass_weighted(Nr_parcels, qr_parcels, rho_water, rho_air, mu, len_parcels)
Nr_array[0,:] = equivalent_number_concentration(D_samples, qr_parcels, rho_water, rho_air)/multiplicity
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

print(sum(Nr_array[0,:]))

# Save with unlimited time dimension
ds.to_netcdf("rain_initiation_stochastic_dataset.nc", unlimited_dims=["time"])
