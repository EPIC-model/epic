import numpy as np
from tools.nc_parcels import nc_parcels
from tools.nc_reader import nc_reader
import subprocess
import os
import argparse

try:
    parser = argparse.ArgumentParser(
            description="Test script to compare serial and parallel merging."
        )

    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=42,
        help="Seeds(s) for the random number generator.",
    )

    parser.add_argument(
        "--n_samples",
        type=int,
        default=42,
        help="Number of random samples per seed.",
    )

    parser.add_argument(
        "--n_ranks",
        type=int,
        nargs='+',
        default=4,
        help="Number of MPI ranks for the parallel runs.",
    )

    parser.add_argument(
        "--n_parcel_per_cell",
        type=int,
        default=30,
        help="Number of parcels (on average) per grid cell.",
    )

    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Print intermediate output."
    )

    mpirun = os.environ.get('MPIRUN')
    exec_path = os.environ.get('EXEC_PATH')

    if mpirun is None:
        raise KeyError('Missing environment variable $MPIRUN.')

    if exec_path is None:
        raise KeyError('Missing environment variable $EXEC_PATH.')


    args = parser.parse_args()

    # nx = ny = nz = 10
    n_parcels = args.n_parcel_per_cell * 10 ** 3

    # domain [0, 1]^3

    #vmin = vcell / parcel%min_vratio
    vcell = 0.001
    min_vratio = 40.0
    vmin = vcell / min_vratio

    tol = 1.0e-14

    exec_parallel = os.path.join(exec_path, 'test_merging_parcels')
    exec_serial = os.path.join(exec_path, 'test_merging_parcels_serial')

    ncp = nc_parcels()
    ncrs = nc_reader()
    ncrp = nc_reader()

    n_fails = 0
    n_merges = 0

    # used when running with --verbose
    total_samples = args.n_samples * len(args.n_ranks) * len(args.seeds)
    modulo = max(1000, int(total_samples / 1000.0))

    i = 0
    for n_rank in args.n_ranks:
        for seed in args.seeds:
            rng = np.random.default_rng(seed)
            for n in range(args.n_samples):
                i = i + 1

                # -------------------------------------------------------------
                # Set up the parcel configuration:
                x_position = rng.uniform(low=0.0, high=1.0, size=n_parcels)
                y_position = rng.uniform(low=0.0, high=1.0, size=n_parcels)
                z_position = rng.uniform(low=0.0, high=1.0, size=n_parcels)

                x_vorticity = rng.uniform(low=-10.0, high=10.0, size=n_parcels)
                y_vorticity = rng.uniform(low=-10.0, high=10.0, size=n_parcels)
                z_vorticity = rng.uniform(low=-10.0, high=10.0, size=n_parcels)

                buoyancy = rng.uniform(low=-1.0, high=1.0, size=n_parcels)

                volume = rng.uniform(low=0.5 * vmin, high=1.5 * vmin, size=n_parcels)

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

                # -------------------------------------------------------------
                # Write the initial NetCDF file:
                ncp.open('initial_parcels.nc')

                ncp.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [10, 10, 10])

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

                # -------------------------------------------------------------
                # Run the serial and parallel versions of the nearest + merging algorithm:
                # We wait 2 minutes per process. If they exceed this time limit, we
                # assume the nearest algorithm is in an endless loop, i.e. deadlocked.
                failed = False
                try:
                    subprocess.run(args=mpirun + ' -np ' + str(n_rank) + ' ' + exec_parallel,
                                   shell=True,
                                   check=True,
                                   timeout=120,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError:
                    n_fails = n_fails + 1
                    failed = True

                try:
                    subprocess.run(args=mpirun + ' -np 1 ' + exec_serial,
                                   shell=True,
                                   check=True,
                                   timeout=120,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError:
                    n_fails = n_fails + 1
                    failed = True

                if not os.path.exists('serial_final_0000000001_parcels.nc') and \
                   not os.path.exists('parallel_final_0000000001_parcels.nc'):
                    failed = True

                if not failed:
                    # ---------------------------------------------------------
                    # Compare the results:
                    ncrs.open('serial_final_0000000001_parcels.nc')
                    ncrp.open('parallel_final_0000000001_parcels.nc')

                    ind1 = np.lexsort((ncrs.get_dataset(step=0, name='x_position'),
                                    ncrs.get_dataset(step=0, name='y_position'),
                                    ncrs.get_dataset(step=0, name='z_position')))

                    ind2 = np.lexsort((ncrp.get_dataset(step=0, name='x_position'),
                                    ncrp.get_dataset(step=0, name='y_position'),
                                    ncrp.get_dataset(step=0, name='z_position')))

                    n_merges = n_merges + n_parcels - ncrs.get_num_parcels(step=0)

                    failed = (not ncrs.get_num_parcels(step=0) == ncrp.get_num_parcels(step=0))

                    if not failed:
                        for attr in ['x_position', 'y_position', 'z_position',
                                    'x_vorticity', 'y_vorticity', 'z_vorticity',
                                    'buoyancy', 'volume', 'B11', 'B12', 'B13', 'B22', 'B23']:
                            ds1 = ncrs.get_dataset(step=0, name=attr)
                            ds2 = ncrp.get_dataset(step=0, name=attr)
                            failed = (failed or max(abs(ds1[ind1] - ds2[ind2])) > tol)

                    ncrs.close()
                    ncrp.close()

                # -------------------------------------------------------------
                # Do clean up:
                if failed:
                    failed = False
                    n_fails = n_fails + 1
                    n_str = str(n_fails).zfill(10)
                    if os.path.exists('initial_parcels.nc'):
                        os.rename('initial_parcels.nc',
                                  'initial_' + n_str + '_parcels.nc')
                    if os.path.exists('serial_final_0000000001_parcels.nc'):
                        os.rename('serial_final_0000000001_parcels.nc',
                                  'serial_fail_' + n_str + '_parcels.nc')
                    if os.path.exists('parallel_final_0000000001_parcels.nc'):
                        os.rename('parallel_final_0000000001_parcels.nc',
                                  'parallel_fail_' + n_str + '_parcels.nc')
                else:
                    os.remove('initial_parcels.nc')
                    os.remove('serial_final_0000000001_parcels.nc')
                    os.remove('parallel_final_0000000001_parcels.nc')

                # -------------------------------------------------------------
                # Intermediate info:
                if args.verbose and (i % modulo == 0):
                    print("#samples, #fails, #merges: ", i, n_fails, n_merges)

    # -------------------------------------------------------------------------
    # Print summary:
    print("--------------------------------------------------------------------")
    print("Total number of samples:      ", args.n_samples * len(args.seeds) * len(args.n_ranks))
    print("Seeds:                        ", args.seeds)
    print("MPI ranks:                    ", args.n_ranks)
    print("Number of parcels per sample: ", n_parcels)
    print("Number of parcels per cell:   ", args.n_parcel_per_cell)
    print("Number of fails:              ", n_fails)
    print("Number of merges:             ", n_merges)

except Exception as ex:
    print(ex)
