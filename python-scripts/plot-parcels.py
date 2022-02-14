#!/usr/bin/env python
import argparse
from tools.plots import plot_parcels
import os
import sys

has_bokeh = True

try:
    from tools.bokeh_plots import bokeh_plot
except:
    has_bokeh = False

try:
    parser = argparse.ArgumentParser(
        description="Plot the parcels of several or individual time steps."
    )

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "--filename", type=str, required=True, help="NetCDF output file of EPIC"
    )

    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        required=False,
        default=0,
        help="steps to plot, " + "e.g. 0 10  100",
    )

    parser.add_argument(
        "--show",
        required=False,
        action="store_true",
        help="show plot instead of saving",
    )

    parser.add_argument(
        "--coloring",
        type=str,
        required=False,
        default="aspect-ratio",
        help="how to color the ellipses",
    )

    parser.add_argument(
        "--fmt",
        type=str,
        required=False,
        default="png",
        help="save format (default: png)",
    )

    parser.add_argument(
        "--xmin",
        type=float,
        required=False,
        default=None,
        help=" float to determine x min",
    )

    parser.add_argument(
        "--xmax",
        type=float,
        required=False,
        default=None,
        help=" float to determine x max",
    )

    parser.add_argument(
        "--ymin",
        type=float,
        required=False,
        default=None,
        help="float to determine y min",
    )

    parser.add_argument(
        "--ymax",
        type=float,
        required=False,
        default=None,
        help="float to determine y max",
    )

    if has_bokeh:
        parser.add_argument(
            "--use-bokeh", required=False, action="store_true", help="use Bokeh to plot"
        )

        parser.add_argument(
            "--display",
            type=str,
            required=False,
            default="full HD",
            help="display (bokeh only)",
        )

        parser.add_argument(
            "--hybrid",
            required=False,
            action="store_true",
            help="hybrid plot (bokeh only)",
        )

    if not "--filename" in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")

    for step in args.steps:
        if has_bokeh and args.use_bokeh:
            bokeh_plot(
                fname=args.filename,
                step=step,
                show=args.show,
                fmt=args.fmt,
                coloring=args.coloring,
                display=args.display,
                xmin=args.xmin,
                xmax=args.xmax,
                ymin=args.ymin,
                ymax=args.ymax,
                hybrid=args.hybrid,
            )
        else:
            plot_parcels(
                fname=args.filename,
                step=step,
                show=args.show,
                fmt=args.fmt,
                coloring=args.coloring,
                xmin=args.xmin,
                xmax=args.xmax,
                ymin=args.ymin,
                ymax=args.ymax,
            )

except Exception as ex:
    print(ex)
