import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tools.h5_reader import H5Reader
import matplotlib.colors as cls
from tools.plots import _plot_ellipses
from tools.plot_beautify import *
import progressbar

class EllipseAnimation:
    """
        22 March 2021
        https://stackoverflow.com/questions/19981054/animating-patch-objects-in-python-matplotlib
        https://stackoverflow.com/questions/44620263/how-to-get-animated-patches-instead-of-n-times-plotted-patches-using-python-and
        https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot
        https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html#matplotlib.animation.FuncAnimation
        https://stackoverflow.com/questions/23791922/saving-animation-with-matplotlib
        https://github.com/niltonvolpato/python-progressbar
    """

    def __init__(self):
        self.ani = None

    def create(self, fname, coloring='aspect-ratio'):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        self.nsteps = self.h5reader.get_num_steps()
        self.extent = self.h5reader.get_mesh_extent()
        self.origin = self.h5reader.get_mesh_origin()
        self.coloring = coloring

        fig = plt.figure(figsize=(12, 4), dpi=300)
        self.ax = plt.gca()

        if coloring == 'aspect-ratio':
            self.vmin = 1.0
            self.vmax = self.h5reader.get_parcel_info('lambda')
        else:
            self.vmin, self.vmax = self.h5reader.get_parcel_min_max(coloring)

        self.norm = cls.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.cmap = plt.cm.viridis_r

        self.ani = animation.FuncAnimation(fig       = fig,
                                           func      = self._update,
                                           init_func = self._init,
                                           frames    = self.nsteps,
                                           interval  = 200,                # Delay between frames in milliseconds
                                           blit      = True,
                                           repeat    = True)


    def save(self, fname, fps=3, dpi=200, extra_args=['-vcodec', 'libx264']):
        self.ani.save(fname, fps=fps, dpi=dpi, extra_args=extra_args)
        print("Save animation in", fname)
        self.h5reader.close()


    def _resize(self):
        # make plot domain 5 percent larger
        self.ax.set_aspect('equal', 'box')
        #self.ax.set_xlim([self.origin[0] - 0.05 * self.extent[0],
                          #self.origin[0] + 1.05 * self.extent[0]])

        #self.ax.set_ylim([self.origin[1] - 0.05 * self.extent[0],
                          #self.origin[0] + 1.05 * self.extent[1]])


    def _update(self, step):

        if (step == 0):
            print("Start processing ...")
            self.bar = progressbar.ProgressBar(maxval=self.nsteps).start()

        self.bar.update(step+1)

        self.ax.clear()
        self._resize()

        _plot_ellipses(self.ax, self.h5reader, step, self.coloring,
                       self.vmin, self.vmax, draw_cbar=(step == 0))

        if step == self.nsteps - 1:
            self.bar.finish()
            print("done.")

        return []

    def _init(self):
        self._resize()
        return []
