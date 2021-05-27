#!/usr/bin/env python
import argparse
from tools.plots import plot_parcels
import os
import sys

try:
    parser = argparse.ArgumentParser(
        description="Plot the parcels of several or individual time steps.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument("--filename",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    parser.add_argument("--steps",
                        type=str,
                        required=False,
                        default='-1',
                        help="steps to plot, a range is specified with a colon, " \
                            + "e.g. 0:10 plots step 0 to 10) (default: -1 [all])")

    parser.add_argument("--show",
                        required=False,
                        action='store_true',
                        help="show plot instead of saving")

    parser.add_argument("--coloring",
                        type=str,
                        required=False,
                        help="how to color the ellipses")


    parser.add_argument("--fmt",
                        type=str,
                        required=False,
                        default="png",
                        help="save format (default: png)")

    if not '--filename' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    steps = args.steps

    if ':' in steps:
        steps = steps.split(':')
        begin=int(steps[0])
        end=int(steps[1])
    elif steps.isdigit():
        end = int(steps)
        begin = end
        if end == -1:
            begin = 0
    else:
        raise IOError("Invalid --steps argument format.")

    if end > -1 and end < begin:
        raise ValueError("Invalid plot range: (begin) " + str(begin) + " > " + str(end) + " (end).")

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")

    plot_parcels(fname=args.filename, begin=begin, end=end, show=args.show,
                  fmt=args.fmt, coloring=args.coloring)

except Exception as ex:
    print(ex)
