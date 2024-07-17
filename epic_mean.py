from pylab import *
import xarray as xr

roots=['fomex','fomex_sfc','fomex_thermo']

for root in roots:
    ds=xr.open_dataset(root+'_fields.nc')
    ds['mfcup']=ds['z_velocity']*(ds['z_velocity']>0.)*(ds['ql']>1e-8)
    ds_mean=ds.mean(dim=["x", "y"])
    ds_mean.to_netcdf(root+'_mean.nc')
    mfmean=ds_mean['mfcup'][-15:-1].mean(dim="t")
    mfmean.to_netcdf(root+'_mfmean.nc')
