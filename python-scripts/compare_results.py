#!/usr/bin/env python
import argparse
from tools.nc_reader import nc_reader
import numpy as np

try:
    parser = argparse.ArgumentParser(
        description="Compare EPIC stats, field and parcel files."
    )

    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "--basename1", type=str, required=True, help="Basepath and basename of first result"
    )

    required.add_argument(
        "--basename2", type=str, required=True, help="Basepath and basename of first result"
    )

    args = parser.parse_args()


    def do_stats_comparison(f1, f2, dsets):

        ncreader1 = nc_reader()
        ncreader1.open(f1)

        ncreader2 = nc_reader()
        ncreader2.open(f2)

        nsteps = ncreader1.get_num_steps()

        if not nsteps == ncreader2.get_num_steps():
            print("Wrong step number")
            exit()

        print(" ... we are checking", nsteps, "steps.")

        max_error = -1 * np.ones((len(dsets)))


        for i in range(nsteps):
            for j, dset in enumerate(dsets):
                d1 = ncreader1.get_dataset(i, dset)
                d2 = ncreader2.get_dataset(i, dset)
                max_error[j] = max(abs(d1 - d2), max_error[j])

        ncreader1.close()
        ncreader2.close()

        for j, dset in enumerate(dsets):
            print(dset.ljust(50), max_error[j])
        return max_error.max()


    def do_dset_comparison(f1, f2, dsets):

        ncreader1 = nc_reader()
        ncreader1.open(f1)

        ncreader2 = nc_reader()
        ncreader2.open(f2)

        nsteps = ncreader1.get_num_steps()

        if not nsteps == ncreader2.get_num_steps():
            print("Wrong step number")
            exit()

        print(" ... we are checking", nsteps, "steps.")

        max_error = -1 * np.ones((len(dsets)))

        for i in range(nsteps):
            for j, dset in enumerate(dsets):
                d1 = ncreader1.get_dataset(i, dset).flatten()
                d2 = ncreader2.get_dataset(i, dset).flatten()
                max_error[j] = max(sum(abs(d1 - d2)), max_error[j])

        ncreader1.close()
        ncreader2.close()

        for j, dset in enumerate(dsets):
            print(dset.ljust(50), max_error[j])
        return max_error.max()


    print("Parcel stats file", end="")

    max_err_res = 0.0

    dsets = ['t',
            'pe',
            'ke',
            'te',
            'n_parcels',
            'n_small_parcel',
            'avg_lam',
            'std_lam',
            'avg_vol',
            'std_vol',
            'x_rms_vorticity',
            'y_rms_vorticity',
            'z_rms_vorticity'
            ]


    tmp = do_stats_comparison(args.basename1 + '_parcel_stats.nc',
                              args.basename2 + '_parcel_stats.nc', dsets)

    max_err_res = max(max_err_res, tmp)

    print()
    print("Field stats file", end="")

    dsets = ['t',
            'rms_v',
            'abserr_v',
            'max_npar',
            'min_npar',
            'avg_npar',
            'avg_nspar'
            ]

    tmp = do_stats_comparison(args.basename1 + '_field_stats.nc',
                              args.basename2 + '_field_stats.nc', dsets)

    max_err_res = max(max_err_res, tmp)

    print()
    print("Field file", end="")

    dsets = ['x',
            'y',
            'z',
            'x_velocity',
            'y_velocity',
            'z_velocity',
            'x_vorticity',
            'y_vorticity',
            'z_vorticity',
            'buoyancy',
            'dry_buoyancy',
            'liquid_water_content'
            ]

    tmp = do_dset_comparison(args.basename1 + '_fields.nc',
                             args.basename2 + '_fields.nc', dsets)

    max_err_res = max(max_err_res, tmp)


    print()
    print("Parcel file", end="")


    dsets = [
            'x_position',
            'y_position',
            'z_position',
            'B11',
            'B12',
            'B13',
            'B22',
            'B23',
            'volume',
            'x_vorticity',
            'y_vorticity',
            'z_vorticity',
            'buoyancy',
            'humidity'
            ]

    # we only need to give one parcel file --> nc_reader will loop over all of them
    tmp = do_dset_comparison(args.basename1 + '_0000000002_parcels.nc',
                             args.basename2 + '_0000000002_parcels.nc', dsets)

    max_err_res = max(max_err_res, tmp)

    print("=" * 54)
    print("Maximum absolute error:".ljust(50), max_err_res)

except Exception as ex:
    print(ex)

