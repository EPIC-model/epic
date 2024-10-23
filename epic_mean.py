from pylab import *
import xarray as xr

ds=xr.open_dataset('fomex_fields.nc')
ds['mfcup']=ds['z_velocity']*(ds['z_velocity']>0.)*(ds['ql']>1e-8)
ds_mean=ds.mean(dim=["x", "y"])
ds_mean.to_netcdf('fomex_mean.nc')
mfmean=ds_mean['mfcup'][-15:-1].mean(dim="t")
mfmean.to_netcdf('fomex_mfmean.nc')
