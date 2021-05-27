from tools.h5_reader import H5Reader
from tools.plot_beautify import *
from tools.plot_style import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os

def _plot_parcels(ax, h5reader, step, coloring, vmin, vmax, draw_cbar=True):

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_mesh_origin()
    extent = h5reader.get_mesh_extent()

    if coloring == 'aspect-ratio':
        data = h5reader.get_aspect_ratio(step=step)
    else:
        data = h5reader.get_parcel_dataset(step=step, name=coloring)

    ells = h5reader.get_ellipses(step=step)
    for j, e in enumerate(ells):
        ax.add_artist(e)
        #e.set_clip_box(ax.bbox)
        e.set_alpha(0.75)
        e.set_facecolor(cmap(norm(data[j])))

    ax.axvline(origin[0], color='black', linewidth=0.25)
    ax.axvline(origin[0] + extent[0], color='black', linewidth=0.25)

    ax.axhline(origin[1], color='black', linewidth=0.25)
    ax.axhline(origin[1] + extent[1], color='black', linewidth=0.25)

    # 26 May 2021
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axis_equal_demo.html
    ax.set_aspect('equal', 'box')

    add_timestamp(ax, h5reader.get_step_attribute(step=step, name='t'))

    add_number_of_parcels(ax, len(data))

    if draw_cbar:
        # 27 May 2021
        # https://stackoverflow.com/questions/29516157/set-equal-aspect-in-plot-with-colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        #fig.add_axes(cax)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, drawedges=False, ax=ax, cax=cax)
        # 19 Feb 2021
        # https://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
        cbar.set_alpha(0.75)
        cbar.solids.set_edgecolor("face")
        cbar.draw_all()

        if coloring == 'aspect ratio':
            cbar.set_label(r'$1 \leq \lambda \leq \lambda_{\max}$')
        else:
            cbar.set_label(coloring)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')


def plot_parcels(fname, begin=0, end=-1, show=False, fmt="png",
                  coloring='aspect-ratio'):
    h5reader = H5Reader()

    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    if end == -1:
        end = nsteps + 1

    if coloring == 'aspect-ratio':
        vmin = 1.0
        vmax = h5reader.get_parcel_info('lambda')
    else:
        vmin, vmax = h5reader.get_parcel_min_max(coloring)

    for i in range(begin, end+1):
        fig, ax = plt.subplots(figsize=(15, 14), dpi=300, num=i)

        _plot_parcels(ax, h5reader, i, coloring, vmin, vmax)

    h5reader.close()

    if show:
        plt.show()
        plt.close()
    else:
        plt.savefig('parcels_'  + coloring + '_step_' + str(i).zfill(len(str(nsteps))) + '.' + fmt,
                    bbox_inches='tight')



def plot_ellipse_orientation(fname, step=0, parcel=0, show=False, fmt="png"):
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


def plot_volume_symmetry_error(fname, show=False, fmt="png"):
    """
    Plot the symmetry error of the gridded volume.
    """
    h5reader = H5Reader()
    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    vmean = np.zeros(nsteps)
    vmin = np.zeros(nsteps)
    vmax = np.zeros(nsteps)
    vstd = np.zeros(nsteps)
    iters = range(nsteps)
    for i in iters:
        volg = h5reader.get_field_dataset(i, 'volume')
        vmean[i] = volg.mean()
        vstd[i] = volg.std()
        vmin[i] = volg.min()
        vmax[i] = volg.max()

    h5reader.close()

    plt.figure()
    plt.grid(which='both', linestyle='dashed')
    #plt.plot(iters, vmean, label='mean', color='darkorange', linewidth=1.5)
    plt.fill_between(iters, vmean - vstd, vmin, label='min-max', color='cornflowerblue',
                     edgecolor='cornflowerblue', linewidth=0.75)
    plt.fill_between(iters, vmean + vstd, vmax, color='cornflowerblue',
                     edgecolor='cornflowerblue', linewidth=0.75)
    plt.fill_between(iters, vmean - vstd, vmean + vstd, label='std', color='royalblue',
                     edgecolor='royalblue', linewidth=0.75)

    plt.xlabel('number of iterations')
    plt.ylabel('volume symmetry error')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_volsymerr.' + fmt, bbox_inches='tight')
    plt.close()


def plot_rms_volume_error(fnames, show=False, fmt="png"):
    """
    Plot the gridded rms volume error.
    """
    n = len(fnames)

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    plt.figure()

    h5reader = H5Reader()
    for i, fname in enumerate(fnames):
        h5reader.open(fname)
        vrms = h5reader.get_diagnostic('rms volume error')
        h5reader.close()

        prefix = os.path.splitext(fname)[0]
        plt.plot(vrms, label=r'' + prefix, linewidth=2, color=colors[i])

    plt.xlabel(r'number of iterations')
    plt.ylabel(r'rms volume error')
    plt.grid(linestyle='dashed', zorder=-1)
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('rms_vol_err.' + fmt, bbox_inches='tight')
    plt.close()


def plot_max_volume_error(fnames, show=False, fmt="png"):
    """
    Plot the gridded absolute volume error (normalised with
    cell volume).
    """
    n = len(fnames)

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    plt.figure()

    h5reader = H5Reader()
    for i, fname in enumerate(fnames):
        h5reader.open(fname)
        vrms = h5reader.get_diagnostic('max absolute normalised volume error')
        h5reader.close()

        prefix = os.path.splitext(fname)[0]
        plt.plot(vrms, label=r'' + prefix, linewidth=2, color=colors[i])

    plt.xlabel(r'number of iterations')
    plt.ylabel(r'max normalised volume error')
    plt.grid(linestyle='dashed', zorder=-1)
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('max_normalised_vol_err.' + fmt, bbox_inches='tight')
    plt.close()


def plot_aspect_ratio(fname, show=False, fmt="png"):
    """
    Plot the mean and standard deviation of the parcel aspect ratio.
    """
    h5reader = H5Reader()
    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    lam_mean = np.zeros(nsteps)
    lam_std = np.zeros(nsteps)

    for step in range(nsteps):
        lam = h5reader.get_aspect_ratio(step)

        lam_mean[step] = lam.mean()
        lam_std[step] = lam.std()

    lmax = h5reader.get_parcel_info('lambda')[0]

    h5reader.close()

    plt.figure()
    plt.plot(lam_mean, color='blue', label=r'mean')
    plt.fill_between(range(nsteps), lam_mean - lam_std, lam_mean + lam_std,
                     alpha=0.5, label=r'std. dev.')

    plt.axhline(lmax, linestyle='dashed', color='black',
                label=r'$\lambda\le\lambda_{\max} = ' + str(lmax) + '$')

    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'number of iterations')
    plt.ylabel(r'aspect ratio $\lambda$')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.25))

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_aspect_ratio_profile.' + fmt, bbox_inches='tight')
    plt.close()


def plot_parcel_number(fname, show=False, fmt="png"):
    """
    Plot the number of parcels in simulation.
    """
    h5reader = H5Reader()
    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    nparcels = np.zeros(nsteps)

    for step in range(nsteps):
        nparcels[step] = h5reader.get_num_parcels(step)

    h5reader.close()

    plt.figure()
    plt.plot(nparcels, color='blue')
    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'number of iterations')
    plt.ylabel(r'parcel count')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_parcel_number_profile.' + fmt, bbox_inches='tight')
    plt.close()


def plot_parcel_volume(fname, show=False, fmt="png"):
    """
    Plot the mean and standard deviation of the parcel volume
    normalised with the cell volume.
    """
    h5reader = H5Reader()
    h5reader.open(fname)

    nsteps = h5reader.get_num_steps()

    vol_mean = np.zeros(nsteps)
    vol_std = np.zeros(nsteps)

    extent = h5reader.get_mesh_extent()
    grid   = h5reader.get_mesh_grid()
    vcell = np.prod(extent / (grid - 1))

    for step in range(nsteps):
        vol = h5reader.get_parcel_dataset(step, 'volume')

        vol_mean[step] = vol.mean() / vcell
        vol_std[step] = vol.std() / vcell

    h5reader.close()

    plt.figure()
    plt.plot(vol_mean, color='blue', label=r'mean')
    plt.fill_between(range(nsteps), vol_mean - vol_std, vol_mean + vol_std,
                     alpha=0.5, label=r'std. dev.')

    plt.axhline(1.0, linestyle='dashed', color='black',
                label=r'cell volume $V_{0}$')

    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'number of iterations')
    plt.ylabel(r'parcel volume / $V_{0}$')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.25))

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_parcel_volume_profile.' + fmt, bbox_inches='tight')
    plt.close()
