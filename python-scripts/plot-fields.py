#!/usr/bin/env python
import argparse
from tools.plots import plot_volume_symmetry_error
import os
import sys

try:
    parser = argparse.ArgumentParser(
        description="Plot the ellipses of several or individual time steps.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    kinds = [
        'volume-symmetry-error'
    ]


    required.add_argument("--filename",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    parser.add_argument("--kind",
                        type=str,
                        required=True,
                        default=kinds[0],
                        help="What kind of field plot: " + str(kinds))

    parser.add_argument("--show",
                        required=False,
                        action='store_true',
                        help="show plot instead of saving")

    parser.add_argument("--fmt",
                        type=str,
                        required=False,
                        default="png",
                        help="save format (default: png)")

    if not '--filename' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")

    if args.kind not in kinds:
        raise ValuError("Plot '" + args.kind + "' not supported!")

    if args.kind == kinds[0]:
        plot_volume_symmetry_error(args.filename, show=args.show, fmt=args.fmt)


except Exception as ex:
    print(ex)
