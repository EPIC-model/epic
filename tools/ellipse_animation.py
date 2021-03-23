#!/usr/bin/env python
import argparse
from animate import EllipseAnimation
import os

try:
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="print some more information")
    parser.add_argument("-f", "--filename",
                        type=str,
                        default='',
                        help="hdf5 file")

    parser.add_argument("-s", "--saveas",
                        type=str,
                        default='',
                        help="file name of saved animation")

    args = parser.parse_args()

    if args.filename == '':
        parser.print_help()
        exit(0)

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")


    if args.saveas == '':
        args.saveas = os.path.split(args.filename)[0] + '.mp4'

    anim = EllipseAnimation(verbose=args.verbose)

    anim.create(args.filename)

    anim.save(args.saveas)

except Exception as ex:
    print(ex)
