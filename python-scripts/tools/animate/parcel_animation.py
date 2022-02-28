import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tools.nc_reader import nc_reader
import matplotlib.colors as cls
from tools.plots import _plot_parcels
from tools.mpl_beautify import *
import progressbar


class ParcelAnimation:
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

    def create(self, fname, coloring="aspect-ratio", **kwargs):
        self.ncreader = nc_reader()

        self.ncreader.open(fname)

        if not self.ncreader.is_parcel_file:
            raise IOError("Not a parcel output file.")

        self.nsteps = self.ncreader.get_num_steps()
        self.coloring = coloring

        self.fkwargs = kwargs

        fig = plt.figure(figsize=(16, 9), dpi=200)
        self.ax = plt.gca()

        if coloring == "aspect-ratio":
            self.vmin = 1.0
            self.vmax = self.ncreader.get_global_attribute("lambda_max")
        else:
            self.vmin, self.vmax = self.ncreader.get_dataset_min_max(coloring)

        self.norm = cls.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.cmap = plt.cm.viridis_r

        self.ani = animation.FuncAnimation(
            fig=fig,
            func=self._update,
            init_func=self._init,
            frames=self.nsteps,
            interval=200,  # Delay between frames in milliseconds
            blit=True,
            repeat=True,
        )

    def save(self, fname, fps=3, dpi=200, extra_args=["-vcodec", "libx264"]):
        self.ani.save(fname, fps=fps, dpi=dpi, extra_args=extra_args)
        print("Save animation in", fname)
        self.ncreader.close()

    def _resize(self):
        # make plot domain 5 percent larger
        self.ax.set_aspect("equal", "box")

    def _update(self, step):

        if step == 0:
            print("Start processing ...")
            self.bar = progressbar.ProgressBar(maxval=self.nsteps).start()

        self.bar.update(step + 1)

        self.ax.clear()
        self._resize()

        _plot_parcels(
            self.ax,
            self.ncreader,
            step,
            self.coloring,
            self.vmin,
            self.vmax,
            draw_cbar=(step == 0),
            **self.fkwargs
        )

        if step == self.nsteps - 1:
            self.bar.finish()
            print("done.")

        return []

    def _init(self):
        self._resize()
        return []
