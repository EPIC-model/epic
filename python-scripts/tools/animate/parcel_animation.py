import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tools.h5_reader import H5Reader
import matplotlib.colors as cls
from tools.plot_beautify import *
import progressbar


class ParcelAnimation:
    """
    14 May 2021
    https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
    """
    def __init__(self):
        self.ani = None

    def create(self, fname):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        self.nsteps = self.h5reader.get_num_steps()
        self.extent = self.h5reader.get_mesh_extent()
        self.origin = self.h5reader.get_mesh_origin()

        fig = plt.figure(figsize=(5, 4), dpi=200)
        self.ax = plt.gca()

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
        self.ax.set_xlim([self.origin[0] - 0.05 * self.extent[0],
                          self.origin[0] + 1.05 * self.extent[0]])

        self.ax.set_ylim([self.origin[1] - 0.05 * self.extent[0],
                          self.origin[0] + 1.05 * self.extent[1]])


    def _update(self, step):

        if (step == 0):
            print("Start processing ...")
            self.bar = progressbar.ProgressBar(maxval=self.nsteps).start()

        self.bar.update(step+1)

        self.ax.clear()

        self._resize()

        self.ax.set_xlabel(r'$x$')
        self.ax.set_ylabel(r'$y$')

        pos = self.h5reader.get_parcel_dataset(step, 'position')

        sc = self.ax.scatter(pos[0, :], pos[1, :], marker='.', s=2)

        self.ax.axvline(self.origin[0], color='black', linewidth=0.25)
        self.ax.axvline(self.origin[0] + self.extent[0], color='black', linewidth=0.25)

        self.ax.axhline(self.origin[1], color='black', linewidth=0.25)
        self.ax.axhline(self.origin[1] + self.extent[1], color='black', linewidth=0.25)

        # add some additional information
        add_timestamp(self.ax, self.h5reader.get_step_attribute(step, name='t'))

        add_number_of_parcels(self.ax, pos.shape[1])

        if step == self.nsteps - 1:
            self.bar.finish()
            print("done.")

        return sc,

    def _init(self):
        self._resize()
        return []
