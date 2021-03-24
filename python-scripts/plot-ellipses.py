#!/usr/bin/env python
import argparse
from tools.plots import plot_ellipses
import os
import sys

try:
    parser = argparse.ArgumentParser(
        description="Plot the ellipses of several or individual time steps.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument("--filename",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    parser.add_argument("--step",
                        type=int,
                        required=False,
                        default=-1,
                        help="step to plot (default: -1 [all])")

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

    plot_ellipses(fname=args.filename, step=args.step, show=args.show, fmt=args.fmt)

except Exception as ex:
    print(ex)
