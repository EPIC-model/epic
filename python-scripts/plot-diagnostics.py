#!/usr/bin/env python
import argparse
from tools.plots import     \
    plot_rms_volume_error,  \
    plot_max_volume_error,  \
    plot_aspect_ratio
import os
import sys

try:
    parser = argparse.ArgumentParser(
        description="Creates diagnostic plots.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    kinds = [
        'rms-volume-error',
        'max-volume-error',
        'aspect-ratio'
    ]


    required.add_argument("--filenames",
                          type=str,
                          nargs='+',
                          required=True,
                          help="list of hdf5 output files of EPIC")

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

    if not '--filenames' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    for fname in args.filenames:
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")

    if args.kind == kinds[0]:
        plot_rms_volume_error(args.filenames, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[1]:
        plot_max_volume_error(args.filenames, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[2]:
        for fname in args.filenames:
            plot_aspect_ratio(fname, show=args.show, fmt=args.fmt)
    else:
        raise ValuError("Plot '" + args.kind + "' not supported!")


except Exception as ex:
    print(ex)
