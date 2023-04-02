import numpy as np
from tools.nc_parcels import nc_parcels
import subprocess
import os

seed = 42

n_runs = 1

n_parcels =

vmin = 1.0

tol = 1.0e-12

exec_parallel = ''
exec_serial = ''


rng = np.random.default_rng(seed)

ncp = nc_parcels()
ncrs = nc_reader()
ncrp = nc_reader()

for i in range(n_runs):

    # -------------------------------------------------------------------------
    # Set up the parcel configuration:
    x_position = rng.uniform(row=0.0, high=1.0, size=n_parcels)
    y_position = rng.uniform(row=0.0, high=1.0, size=n_parcels)
    z_position = rng.uniform(row=0.0, high=1.0, size=n_parcels)

    x_vorticity = rng.uniform(row=-10.0, high=10.0, size=n_parcels)
    y_vorticity = rng.uniform(row=-10.0, high=10.0, size=n_parcels)
    z_vorticity = rng.uniform(row=-10.0, high=10.0, size=n_parcels)

    buoyancy = rng.uniform(low=-1.0, high=1.0, size=n_parcels)

    volume = rng.uniform(low=0.5 * vmin, high=2.0 * vmin, size=n_parcels)

    # lam = a / c
    lam = rng.uniform(low=1.0, high=4.0, size=n_parcels)

    # lam2 = a / b
    lam2 = rng.uniform(low=1.0, high=4.0, size=n_parcels)

    abc = 0.75 * volume / np.pi

    a2 = (abc * lam * lam2)  ** (2.0/3.0)
    b2 = a2 / lam2 ** 2
    c2 = a2 / lam ** 2

    theta = rng.uniform(low=0.0, high=2.0*np.pi, size=n_parcels)
    phi = rng.uniform(low=0.0, high=2.0*np.pi, size=n_parcels)

    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
    B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
    B13 = (a2 - c2) * ct * sp * cp
    B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
    B23 = (a2 - c2) * st * sp * cp

    # -------------------------------------------------------------------------
    # Write the initial NetCDF file:
    ncp.open('initial_configuration.nc')

    ncp.add_box(tuple_origin, tuple_extent, tuple_ncells)

    ncp.add_dataset('x_position', x_position, unit='m')
    ncp.add_dataset('y_position', y_position, unit='m')
    ncp.add_dataset('z_position', z_position, unit='m')

    ncp.add_dataset('buoyancy', buoyancy, unit='m/s^2')

    ncp.add_dataset('volume', volume, unit='m^3')

    ncp.add_dataset('x_vorticity', x_vorticity, unit='1/s')
    ncp.add_dataset('y_vorticity', y_vorticity, unit='1/s')
    ncp.add_dataset('z_vorticity', z_vorticity, unit='1/s')

    ncp.add_dataset('B11', B11, unit='m^2')
    ncp.add_dataset('B12', B12, unit='m^2')
    ncp.add_dataset('B13', B13, unit='m^2')
    ncp.add_dataset('B22', B22, unit='m^2')
    ncp.add_dataset('B23', B23, unit='m^2')

    ncp.close()

    # -------------------------------------------------------------------------
    # Run the serial and parallel versions of the nearest + merging algorithm:
    process_1 = subprocess.Popen(['mpirun', '-np', n_ranks, exec_parallel])
    process_2 = subprocess.Popen(['mpirun', '-np', 1, exec_serial])

    process_1.wait()
    process_2.wait()

    # -------------------------------------------------------------------------
    # Compare the results:
    ncrs.open('serial_final_configuration.nc')
    ncrp.open('parallel_final_configuration.nc')

    failed = False
    for attr in ['x_position', 'y_position', 'z_position',
                 'x_vorticity', 'y_vorticity', 'z_vorticity',
                 'buoyancy', 'volume', 'B11', 'B12', 'B13', 'B22', 'B23']:
        ds1 = ncrs.get_dataset(step=0, name=attr)
        ds2 = ncrp.get_dataset(step=0, name=attr)
        failed = (failed .or. max(abs(ds1 - ds2)) > tol)

    ncrs.close()
    ncrp.close()

    if failed:
        os.rename('initial_configuration.nc', 'initial_configuration_' + str(i) + '.nc')
        os.rename('serial_final_configuration.nc', 'serial_final_configuration_' + str(i) + '.nc')
        os.rename('parallel_final_configuration.nc', 'parallel_final_configuration_' + str(i) + '.nc')
        failed = False
