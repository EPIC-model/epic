#!/usr/bin/env python
import argparse
from tools.animate import EllipseAnimation
import os

try:
    parser = argparse.ArgumentParser(
        description="Save a mp4 animation of the evolving ellipses.")

    parser.add_argument("-f", "--filename",
                        type=str,
                        default='',
                        help="hdf5 output file of EPIC")

    parser.add_argument("-s", "--saveas",
                        type=str,
                        default='',
                        help="file name of saved animation (default: FILENAME.mp4)")

    args = parser.parse_args()

    if args.filename == '':
        parser.print_help()
        exit(0)

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")


    if args.saveas == '':
        args.saveas = os.path.splitext(args.filename)[0] + '.mp4'
    else:
        args.saveas = os.path.splitext(args.saveas)[0] + '.mp4'

    anim = EllipseAnimation()

    anim.create(args.filename)

    anim.save(args.saveas)

except Exception as ex:
    print(ex)
