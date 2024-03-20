import numpy as np
from tools.nc_reader import nc_reader
import subprocess
import os
import argparse

try:
    parser = argparse.ArgumentParser(
            description="Test script to compare serial and parallel merging."
        )

    parser.add_argument(
        "--n_samples",
        type=int,
        default=100,
        help="Number of random samples.",
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
        "--min_vratio",
        type=float,
        default=40.0,
        help="min volume ratio",
    )

    parser.add_argument(
        "--nx",
        type=int,
        default=32,
        help="number of grid cells in x",
    )

    parser.add_argument(
        "--ny",
        type=int,
        default=32,
        help="number of grid cells in y",
    )

    parser.add_argument(
        "--nz",
        type=int,
        default=32,
        help="number of grid cells in z",
    )

    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Print intermediate output."
    )

    parser.add_argument(
        "--cmd",
        type=str,
        default="mpirun",
        help="Run jobs with 'mpirun' or 'srun'."
    )

    exec_path = os.environ.get('EXEC_PATH')

    if exec_path is None:
        raise KeyError('Missing environment variable $EXEC_PATH.')


    args = parser.parse_args()

    if not args.cmd == 'mpirun' and not args.cmd == 'srun':
        raise IOError("Use either 'mpirun' or 'srun'.")

    nx = args.nx
    ny = args.ny
    nz = args.nz
    n_parcels = args.n_parcel_per_cell * nx * ny * nz

    # domain [0, 1]^3

    #vmin = vcell / parcel%min_vratio
    vcell = 1.0 / (nx * ny * nz)
    min_vratio = args.min_vratio
    vmin = vcell / min_vratio

    tol = 1.0e-14

    exe = os.path.join(exec_path, 'verify_parcel_merging')

    flags = ' --nx ' + str(nx) \
          + ' --ny ' + str(ny) \
          + ' --nz ' + str(nz) \
          + ' --min_vratio ' + str(min_vratio)

    ncrs = nc_reader()
    ncrp = nc_reader()

    n_fails = 0
    n_merges = 0

    # used when running with --verbose
    modulo = 100

    ntasks_per_node = 128

    for n in range(args.n_samples):
        try:
            cmd = 'mpirun -np 1 '
            if args.cmd == 'srun':
                cmd = 'srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --exact '
            subprocess.run(args=cmd + exe + flags + ' --setup-parcels',
                           shell=True,
                           check=True,
                           timeout=360,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            print('Error in running the serial version.', flush=True)

        # Note: Because the 1 MPI run writes the initial parcel setup; the actual
        # solve has number 2 instead of 1.
        if os.path.exists('serial_final_0000000002_parcels.nc'):
            os.system('mv serial_final_0000000002_parcels.nc serial_final_0000000001_parcels.nc')
            ncrs.open('serial_final_0000000001_parcels.nc')
            n_merges = n_merges + n_parcels - ncrs.get_num_parcels(step=0)
            ncrs.close()
            print ("Sample", n, "generated.", flush=True)

        for n_rank in args.n_ranks:
            #print("Run sample with", n_rank, "MPI ranks", flush=True)

            nodes = int(np.ceil(n_rank / ntasks_per_node))

            # -------------------------------------------------------------
            # Run the serial and parallel versions of the nearest + merging algorithm:
            # We wait 2 minutes per process. If they exceed this time limit, we
            # assume the nearest algorithm is in an endless loop, i.e. deadlocked.
            failed = False
            try:
                cmd = 'mpirun -np ' + str(n_rank) + ' '
                if args.cmd == 'srun':
                    cmd = 'srun --nodes=' + str(nodes) + ' --ntasks=' + str(n_rank)
                    cmd = cmd + ' --cpus-per-task=1 --exact '
                subprocess.run(args=cmd + exe + flags,
                               shell=True,
                               check=True,
                               timeout=360,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                print('Error in running the parallel version.', flush=True)
                failed = True

            if not os.path.exists('serial_final_0000000001_parcels.nc') and \
                not os.path.exists('parallel_final_0000000001_parcels.nc'):
                failed = True

            if not failed:
                # ---------------------------------------------------------
                # Compare the results:
                ncrs.open('serial_final_0000000001_parcels.nc')
                ncrp.open('parallel_final_0000000001_parcels.nc')

                if ncfile.dimensions['world%size'].size != n_rank:
                    print('Warning: Failed to run with requested number of MPI cores', flush=True)

                ind1 = np.lexsort((ncrs.get_dataset(step=0, name='x_position'),
                                ncrs.get_dataset(step=0, name='y_position'),
                                ncrs.get_dataset(step=0, name='z_position')))

                ind2 = np.lexsort((ncrp.get_dataset(step=0, name='x_position'),
                                ncrp.get_dataset(step=0, name='y_position'),
                                ncrp.get_dataset(step=0, name='z_position')))

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
                if os.path.exists('initial_0000000001_parcels.nc'):
                    os.rename('initial_0000000001_parcels.nc',
                              'initial_' + n_str + '_parcels.nc')
                if os.path.exists('serial_final_0000000001_parcels.nc'):
                    os.system('cp serial_final_0000000001_parcels.nc serial_fail_' + n_str + '_parcels.nc')
                if os.path.exists('parallel_final_0000000001_parcels.nc'):
                    os.rename('parallel_final_0000000001_parcels.nc',
                              'parallel_fail_' + n_str + '_parcels.nc')
            else:
                os.remove('parallel_final_0000000001_parcels.nc')

        # -------------------------------------------------------------
        # Intermediate info:
        if args.verbose and (n % modulo == 0):
            print("#samples, #fails, #merges: ", n, n_fails, n_merges, flush=True)

        os.remove('initial_0000000001_parcels.nc')
        if os.path.exists('serial_final_0000000001_parcels.nc'):
            os.remove('serial_final_0000000001_parcels.nc')

    # -------------------------------------------------------------------------
    # Print summary:
    print("--------------------------------------------------------------------", flush=True)
    print("Total number of samples:      ", args.n_samples, flush=True)
    print("MPI ranks:                    ", args.n_ranks, flush=True)
    print("Number of parcels per sample: ", n_parcels, flush=True)
    print("Number of parcels per cell:   ", args.n_parcel_per_cell, flush=True)
    print("Number of fails:              ", n_fails, flush=True)
    print("Number of merges:             ", n_merges, flush=True)

except Exception as ex:
    print(ex, flush=True)
