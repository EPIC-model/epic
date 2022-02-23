#!/usr/bin/env python
from tools.nc_reader import nc_reader
from tools.mpl_beautify import *
from tools.mpl_style import *
import argparse
import numpy as np
import matplotlib.pyplot as plt

try:

    class IndexTracker:
        """
        Reference: (22 Feb 2022)
            https://matplotlib.org/stable/gallery/event_handling/image_slices_viewer.html
        """
        def __init__(self, ax, ncreader, step, name):
            self.ax = ax
            self.ann = None
            self.bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
            #ax.set_title('use scroll wheel to navigate images')

            # 22 Feb 2022
            # https://stackoverflow.com/questions/32034237/how-does-numpys-transpose-method-permute-the-axes-of-an-array
            self.X = ncreader.get_dataset(step, name).transpose(0, 1, 2)
            rows, cols, self.slices = self.X.shape
            self.ind = self.slices//2

            long_name = ncreader.get_dataset_attribute(name, 'long_name')
            units = ncreader.get_dataset_attribute(name, 'units')

            origin = ncreader.get_box_origin()
            extent = ncreader.get_box_extent()

            self.im = ax.imshow(self.X[:, :, self.ind],
                                origin='lower',
                                extent=(origin[0], origin[0]+extent[0],
                                        origin[2], origin[2]+extent[2]),
                                interpolation='bilinear')

            add_timestamp(self.ax, ncreader.get_dataset(step=step, name="t"),
                          xy=(0.05, 1.05), fmt="%.3f")

            ax.set_xlabel(r'horizontal width ($\si{m}$)')
            ax.set_ylabel(r'height ($\si{m}$)')


            cbar = self.im.axes.figure.colorbar(self.im, ax=ax)
            cbar.set_label(long_name + r' ($\si{' + units + '}$)')
            #self.im.axes.figure.tight_layout()
            self.update()

        def on_scroll(self, event):
            print("%s %s" % (event.button, event.step))
            if event.button == 'up':
                self.ind = (self.ind + 1) % self.slices
            else:
                self.ind = (self.ind - 1) % self.slices
            self.update()

        def update(self):
            self.im.set_data(self.X[:, :, self.ind])
            text = 'slice %s' % self.ind
            # 22 Feb 2022
            # https://stackoverflow.com/questions/42314525/remove-annotation-while-keeping-plot-matplotlib
            if self.ann is not None:
                self.ann.remove()
            self.ann = self.ax.annotate(text, xy=(0.75, 1.05), xycoords="axes fraction", bbox=self.bbox)
            #self.im.axes.figure.tight_layout()
            self.im.axes.figure.canvas.draw()



    parser = argparse.ArgumentParser(
        description="Show slice plots of EPIC fields at given time steps."
    )

    parser.add_argument(
        "--filename",
        type=str,
        required=True,
        help="NetCDF output field file of EPIC"
    )

    parser.add_argument(
        "--step",
        type=int,
        required=True,
        help="step to show",
    )

    parser.add_argument(
        "--field",
        type=str,
        required=True,
        help="field to show",
    )

    args = parser.parse_args()

    fig, ax = plt.subplots(1, 1)

    ncreader = nc_reader()

    ncreader.open(args.filename)
    nsteps = ncreader.get_num_steps()

    tracker = IndexTracker(ax, ncreader, args.step, args.field)

    ncreader.close()

    fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
    plt.show()
except Exception as ex:
    print(ex)
