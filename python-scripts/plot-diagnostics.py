#!/usr/bin/env python
import argparse
from tools.plots import     \
    plot_rms_volume_error,  \
    plot_max_volume_error,  \
    plot_aspect_ratio,      \
    plot_parcel_volume,     \
    plot_parcel_number,     \
    plot_center_of_mass,    \
    plot_cumulative,        \
    plot_parcels_per_cell
import os
import sys

try:
    parser = argparse.ArgumentParser(description="Creates diagnostic plots.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    kinds = [
        'rms-volume-error',
        'max-volume-error',
        'aspect-ratio',
        'parcel-volume',
        'parcel-number',
        'center-of-mass',
        'parcel-cumulative',
        'number-parcel-per-cell'
    ]


    required.add_argument("--filenames",
                          type=str,
                          nargs='+',
                          required=True,
                          help="list of hdf5 output files of EPIC")

    parser.add_argument("--labels",
                        type=str,
                        nargs='+',
                        required=False,
                        help="special labels for the files")

    parser.add_argument("--kind",
                        type=str,
                        required=True,
                        default=kinds[0],
                        help="what kind of diagnostic plot: " + str(kinds))

    parser.add_argument("--show",
                        required=False,
                        action='store_true',
                        help="show plot instead of saving")

    parser.add_argument("--fmt",
                        type=str,
                        required=False,
                        default="png",
                        help="save format (default: png)")

    parser.add_argument("--step",
                        type=int,
                        required=False,
                        default=0,
                        help="step in hdf5 file (cumulative plot only)")


    parser.add_argument("--dataset",
                        type=str,
                        required=False,
                        default='volume',
                        help="parcel attribute (cumulative plot only)")


    if not '--filenames' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    for fname in args.filenames:
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")

    if args.kind == kinds[0]:
        plot_rms_volume_error(args.filenames, show=args.show, fmt=args.fmt,
                              labels=args.labels)
    elif args.kind == kinds[1]:
        plot_max_volume_error(args.filenames, show=args.show, fmt=args.fmt,
                              labels=args.labels)
    elif args.kind == kinds[2]:
        for fname in args.filenames:
            plot_aspect_ratio(fname, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[3]:
        for fname in args.filenames:
            plot_parcel_volume(fname, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[4]:
        plot_parcel_number(args.filenames, show=args.show, fmt=args.fmt,
                           labels=args.labels)
    elif args.kind == kinds[5]:
        for fname in args.filenames:
            plot_center_of_mass(fname, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[6]:
        plot_cumulative(args.filenames, step=args.step,
                        dset=args.dataset, show=args.show,
                        fmt=args.fmt, labels=args.labels)
    elif args.kind == kinds[7]:
        plot_parcels_per_cell(args.filenames,
                              show=args.show,
                              fmt=args.fmt,
                              labels=args.labels)
    else:
        raise ValueError("Plot '" + args.kind + "' not supported!")


except Exception as ex:
    print(ex)
