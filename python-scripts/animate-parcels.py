#!/usr/bin/env python
import argparse
from tools.animate import ParcelAnimation
import os
import sys

has_bokeh = True
try:
    from tools.animate.bokeh_animation import BokehAnimation
except:
    has_bokeh = False

try:
    parser = argparse.ArgumentParser(
        description="Save a mp4 animation of the evolving parcels.")

    # 24 March 2021
    # https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
    required = parser.add_argument_group('required arguments')

    required.add_argument("--filename",
                          type=str,
                          required=True,
                          help="hdf5 output file of EPIC")

    parser.add_argument("-s", "--saveas",
                        type=str,
                        required=False,
                        default='',
                        help="file name of saved animation (default: FILENAME.mp4)")

    parser.add_argument("--coloring",
                        type=str,
                        required=False,
                        default='aspect-ratio',
                        help="how to color the parcels")

    parser.add_argument("--xmin",
                    type=float,
                    required=False,
                    default=None,
                    help=" float to determine x min")

    parser.add_argument("--xmax",
                        type=float,
                        required=False,
                        default=None,
                        help=" float to determine x max")

    parser.add_argument("--ymin",
                        type=float,
                        required=False,
                        default=None,
                        help="float to determine y min")

    parser.add_argument("--ymax",
                        type=float,
                        required=False,
                        default=None,
                        help="float to determine y max")

    if has_bokeh:
        parser.add_argument("--use-bokeh",
                            required=False,
                            action='store_true',
                            help="use Bokeh to plot")


    if not '--filename' in sys.argv:
        parser.print_help()
        exit(0)

    args = parser.parse_args()

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")


    if args.saveas == '':
        args.saveas = os.path.splitext(args.filename)[0] + '-' + args.coloring + '.mp4'
    else:
        args.saveas = os.path.splitext(args.saveas)[0] + '.mp4'

    if has_bokeh and args.use_bokeh:
        anim = BokehAnimation()
    else:
        anim = ParcelAnimation()

    anim.create(args.filename,
                coloring=args.coloring,
                xmin=args.xmin,
                xmax=args.xmax,
                ymin=args.ymin,
                ymax=args.ymax)

    anim.save(args.saveas)

except Exception as ex:
    print(ex)
