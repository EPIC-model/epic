#!/usr/bin/env python
import argparse
from tools.plots import plot_time_pie_chart, plot_time_bar, plot_time_speedup
import os
import sys

try:
    parser = argparse.ArgumentParser(description="Creates timing plots.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group("required arguments")

    kinds = ["pie-chart", "bar-plot", "speedup"]

    required.add_argument(
        "--filenames",
        type=str,
        nargs="+",
        required=True,
        help="list of hdf5 output files of EPIC",
    )

    parser.add_argument(
        "--labels",
        type=str,
        nargs="+",
        required=False,
        help="special labels for the files (bar plot only)",
    )

    parser.add_argument(
        "--nthreads",
        type=int,
        nargs="+",
        required=False,
        help="number of OpenMP threads (speedup plot only)",
    )

    parser.add_argument(
        "--kind",
        type=str,
        required=True,
        default=kinds[0],
        help="what kind of time plot: " + str(kinds),
    )

    parser.add_argument(
        "--show",
        required=False,
        action="store_true",
        help="show plot instead of saving",
    )

    parser.add_argument(
        "--fmt",
        type=str,
        required=False,
        default="png",
        help="save format (default: png)",
    )

    if not "--filenames" in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    for fname in args.filenames:
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")

    if args.kind == kinds[0]:
        for fname in args.filenames:
            plot_time_pie_chart(fname, show=args.show, fmt=args.fmt)
    elif args.kind == kinds[1]:
        plot_time_bar(args.filenames, show=args.show, fmt=args.fmt, labels=args.labels)
    elif args.kind == kinds[2]:
        plot_time_speedup(
            args.filenames, nthreads=args.nthreads, show=args.show, fmt=args.fmt
        )
    else:
        raise ValueError("Plot '" + args.kind + "' not supported!")


except Exception as ex:
    print(ex)
