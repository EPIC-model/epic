from h5_reader import H5Reader
from plot_beautify import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def plot_ellipses(fname, step=-1, show=False, fmt="png"):
    h5reader = H5Reader()

    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    start = 0
    stop = nsteps

    if step > -1:
        start = step
        stop = step + 1

    lam = h5reader.get_parcel_info('lambda')

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=1.0, vmax=lam)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_mesh_origin()
    extent = h5reader.get_mesh_extent()

    for i in range(start, stop):
        ells = h5reader.get_ellipses(step=i)

        fig, ax = plt.subplots(figsize=(5, 4), dpi=300, num=i)

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
            plt.savefig('ellipses_step_' + str(i).zfill(len(str(nsteps))) + '.' + fmt,
                        bbox_inches='tight')
        plt.close()

    h5reader.close()


def plot_ellipse_orientation(fname, step=0, parcel=0, show=False, fmt="png"):
    mpl.rcParams['text.usetex'] = True

    h5reader = H5Reader()

    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    if step > nsteps - 1:
        raise ValueError('Step ' + str(step) + ' > ' + str(nsteps - 1))


    num_parcels = h5reader.get_num_parcels(step)
    if parcel > num_parcels - 1:
        raise ValueError('Parcel index ' + str(parcel) + ' > ' + str(num_parcels - 1))

    lam = h5reader.get_parcel_info('lambda')

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=1.0, vmax=lam)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_mesh_origin()
    extent = h5reader.get_mesh_extent()

    ell = h5reader.get_ellipses(step=step)[parcel]

    vol = h5reader.get_parcel_dataset(step=step, name='volume')[parcel]
    ratio = h5reader.get_aspect_ratio(step=step)[parcel]
    ab = vol / np.pi
    a = 1.02 * np.sqrt(ab * ratio)

    fig, ax = plt.subplots(figsize=(5, 4), dpi=300, num=step)
    ax.axis('off')
    ax.set_aspect('equal')

    angles = np.arange(0, 359, 15)
    radii = np.linspace(0.0, a, 10)
    r = 1.02 * radii[-1]

    # 3 April 2021
    # https://stackoverflow.com/questions/37246941/specifying-the-order-of-matplotlib-layers
    for angle in angles:
        rad = np.deg2rad(angle)
        ax.plot([0, r * np.cos(rad)], [0, r * np.sin(rad)], color='gray',
                linewidth=0.5, zorder=0)
        # 5 April 2021
        # https://stackoverflow.com/questions/19926246/inserting-a-degree-symbol-into-python-plot
        ax.text(1.1 * r * np.cos(rad), 1.1 * r * np.sin(rad), str(angle) + r'$^\circ$',
                horizontalalignment='center', verticalalignment='center')

    angles = np.linspace(0, 2.0 * np.pi, 360)
    for r in radii:
        ax.plot(r * np.cos(angles), r * np.sin(angles), color='gray', linewidth=0.5, zorder=0)

    ax.add_artist(ell)
    ell.set_clip_box(ax.bbox)
    ell.set_alpha(1.0)
    ell.set_facecolor(cmap(norm(ratio)))

    # always use center (0, 0)
    ell.set_center(xy=(0, 0))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, drawedges=False)
    ## 19 Feb 2021
    ## https://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
    cbar.set_alpha(1.0)
    cbar.solids.set_edgecolor("face")
    cbar.draw_all()
    cbar.set_label(r'$1 \leq \lambda \leq \lambda_{\max}$')

    phi = np.rad2deg(h5reader.get_parcel_dataset(step, 'orientation')[parcel])
    add_box(plt, label=r'ellipse orientation', value=phi, unit=r'$^\circ$', fmt='%-3.1f')

    plt.xlim([-1.2 * a, 1.2 * a])
    plt.ylim([-1.2 * a, 1.2 * a])

    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('ellipse_' + str(parcel).zfill(len(str(num_parcels))) + \
                    '_orientation_step_' + str(step).zfill(len(str(nsteps))) + '.' + fmt,
                    bbox_inches='tight')
    plt.close()

    h5reader.close()
