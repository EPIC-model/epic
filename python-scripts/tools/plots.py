from h5_reader import H5Reader
from plot_beautify import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_ellipses(fname, step=-1, show=False):
    h5reader = H5Reader()

    h5reader.open(fname)

    start = 0
    stop = h5reader.get_num_steps()

    if step > -1:
        start = step
        stop = step + 1

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=1.0, vmax=5.0)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_mesh_origin()
    extent = h5reader.get_mesh_extent()

    for i in range(start, stop):
        ells = h5reader.get_ellipses(step=i)

        fig, ax = plt.subplots(figsize=(5, 4), dpi=300)

        ratio = h5reader.get_aspect_ratio(step=i)
        for j, e in enumerate(ells):
            ax.add_artist(e)
            #e.set_clip_box(ax.bbox)
            e.set_alpha(0.75)
            e.set_facecolor(cmap(norm(ratio[j])))

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, drawedges=False)
        ## 19 Feb 2021
        ## https://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
        cbar.set_alpha(0.75)
        cbar.solids.set_edgecolor("face")
        cbar.draw_all()
        cbar.set_label(r'$1 \leq \lambda \leq \lambda_{\max}$')


        plt.axvline(origin[0], color='black', linewidth=0.25)
        plt.axvline(origin[0] + extent[0], color='black', linewidth=0.25)

        plt.axhline(origin[1], color='black', linewidth=0.25)
        plt.axhline(origin[1] + extent[1], color='black', linewidth=0.25)

        add_timestamp(plt, h5reader.get_step_attribute(step=i, name='t'))

        add_number_of_parcels(plt, len(ratio))

        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')
        plt.tight_layout()

        if show:
            plt.show()
        else:
            plt.savefig('ellipses_step_' + str(i).zfill(10) + '.png', bbox_inches='tight')
        plt.close()

    h5reader.close()
