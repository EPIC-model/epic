import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tools.h5_reader import H5Reader
import matplotlib.colors as cls
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

    def create(self, fname):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        self.nsteps = self.h5reader.get_num_steps()

        fig = plt.figure(figsize=(5, 4), dpi=200)
        self.ax = plt.gca()
        self.norm = cls.Normalize(vmin=1.0, vmax=5.0)
        self.cmap = plt.cm.viridis_r

        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.norm)
        cbar = fig.colorbar(sm, drawedges=False)
        # 19 Feb 2021
        # https://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
        cbar.set_alpha(0.75)
        cbar.solids.set_edgecolor("face")
        cbar.draw_all()
        cbar.set_label(r'$1 \leq \lambda \leq \lambda_{\max}$')

        self.ani = animation.FuncAnimation(fig       = fig,
                                           func      = self._update,
                                           init_func = self._init,
                                           frames    = self.nsteps,
                                           interval  = 200,                # Delay between frames in milliseconds
                                           blit      = True)


    def save(self, fname, fps=3, dpi=200, extra_args=['-vcodec', 'libx264']):
        self.ani.save(fname, fps=fps, dpi=dpi, extra_args=extra_args)
        print("Save animation in", fname)
        self.h5reader.close()


    def _update(self, step):

        if (step == 0):
            print("Start processing ...")
            self.bar = progressbar.ProgressBar(maxval=self.nsteps).start()

        self.bar.update(step+1)

        ells = self.h5reader.get_ellipses(step)
        patches = []
        self.ax.clear()
        self.ax.set_xlim([-0.5 * np.pi, 0.5 * np.pi])
        self.ax.set_ylim([-0.5 * np.pi, 0.5 * np.pi])
        self.ax.set_xlabel(r'$x$')
        self.ax.set_ylabel(r'$y$')

        ratio = self.h5reader.get_aspect_ratio(step)
        for j, e in enumerate(ells):
            e.set_alpha(0.75)
            e.set_facecolor(self.cmap(self.norm(ratio[j])))
            patches.append(self.ax.add_patch(e))

        extent = self.h5reader.get_mesh_extent()
        origin = self.h5reader.get_mesh_origin()

        self.ax.axvline(origin[0], color='black', linewidth=0.25)
        self.ax.axvline(origin[0] + extent[0], color='black', linewidth=0.25)

        self.ax.axhline(origin[1], color='black', linewidth=0.25)
        self.ax.axhline(origin[1] + extent[1], color='black', linewidth=0.25)

        add_timestamp(self.ax, self.h5reader.get_step_attribute(step, name='t'))

        add_number_of_parcels(self.ax, len(ratio))
        plt.tight_layout()

        if step == self.nsteps - 1:
            self.bar.finish()
            print("done.")

        return patches

    def _init(self):
        return []
