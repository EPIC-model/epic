#!/usr/bin/env python
import argparse
from tools.animate.bokeh_animation import BokehAnimation
import os
import sys

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

    anim = BokehAnimation()

    anim.create(args.filename, coloring=args.coloring)

    anim.save(args.saveas)

except Exception as ex:
    print(ex)
