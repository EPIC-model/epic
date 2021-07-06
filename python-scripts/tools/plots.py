from tools.h5_reader import H5Reader
from tools.plot_beautify import *
from tools.plot_style import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns

def _plot_parcels(ax, h5reader, step, coloring, vmin, vmax, draw_cbar=True):

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_box_origin()
    extent = h5reader.get_box_extent()

    if coloring == 'aspect-ratio':
        data = h5reader.get_aspect_ratio(step=step)
    else:
        data = h5reader.get_dataset(step=step, name=coloring)

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

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    if end == -1:
        end = nsteps + 1

    if coloring == 'aspect-ratio':
        vmin = 1.0
        vmax = h5reader.get_parcel_option('lambda')
    else:
        vmin, vmax = h5reader.get_dataset_min_max(coloring)

    for i in range(begin, end+1):
        fig, ax = plt.subplots(figsize=(15, 14), dpi=300, num=i)

        _plot_parcels(ax, h5reader, i, coloring, vmin, vmax)

        if show:
            plt.show()
        else:
            plt.savefig('parcels_'  + coloring + '_step_' + str(i).zfill(len(str(nsteps))) + '.' + fmt,
                        bbox_inches='tight')
        plt.close()
    h5reader.close()



def plot_ellipse_orientation(fname, step=0, parcel=0, show=False, fmt="png"):
    h5reader = H5Reader()

    h5reader.open(fname)

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    if step > nsteps - 1:
        raise ValueError('Step ' + str(step) + ' > ' + str(nsteps - 1))


    num_parcels = h5reader.get_num_parcels(step)
    if parcel > num_parcels - 1:
        raise ValueError('Parcel index ' + str(parcel) + ' > ' + str(num_parcels - 1))

    lam = h5reader.get_parcel_option('lambda')

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=1.0, vmax=lam)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_box_origin()
    extent = h5reader.get_box_extent()

    ell = h5reader.get_ellipses(step=step)[parcel]

    vol = h5reader.get_dataset(step=step, name='volume')[parcel]
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

    phi = np.rad2deg(h5reader.get_dataset(step, 'orientation')[parcel])
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
    (The gridded symmetry volume is only written in debug mode.)
    """
    h5reader = H5Reader()
    h5reader.open(fname)

    if h5reader.is_parcel_file:
        raise IOError('Not a field output file.')

    try:
        h5reader.get_dataset(0, 'symmetry volume')
    except:
        raise IOError('This plot is only available in debug mode.')

    nsteps = h5reader.get_num_steps()

    vmean = np.zeros(nsteps)
    vmin = np.zeros(nsteps)
    vmax = np.zeros(nsteps)
    vstd = np.zeros(nsteps)
    t = np.zeros(nsteps)
    for i in range(nsteps):
        sym_volg = h5reader.get_dataset(i, 'symmetry volume')
        vmean[i] = sym_volg.mean()
        vstd[i] = sym_volg.std()
        vmin[i] = sym_volg.min()
        vmax[i] = sym_volg.max()
        t[i] = h5reader.get_step_attribute(i, 't')

    h5reader.close()

    plt.figure()
    plt.grid(which='both', linestyle='dashed')
    #plt.plot(iters, vmean, label='mean', color='darkorange', linewidth=1.5)
    plt.fill_between(t, vmean - vstd, vmin, label='min-max', color='cornflowerblue',
                     edgecolor='cornflowerblue', linewidth=0.75)
    plt.fill_between(t, vmean + vstd, vmax, color='cornflowerblue',
                     edgecolor='cornflowerblue', linewidth=0.75)
    plt.fill_between(t, vmean - vstd, vmean + vstd, label='std', color='royalblue',
                     edgecolor='royalblue', linewidth=0.75)

    plt.xlabel(r'time (s)')
    plt.ylabel(r'volume symmetry error')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_volsymerr.' + fmt, bbox_inches='tight')
    plt.close()


def plot_rms_volume_error(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the gridded rms volume error.
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    plt.figure()

    h5reader = H5Reader()
    for i, fname in enumerate(fnames):
        h5reader.open(fname)

        if h5reader.is_parcel_file:
            raise IOError('Not a field output file.')

        vrms = h5reader.get_diagnostic('rms volume error')
        t = h5reader.get_diagnostic('t')
        h5reader.close()

        label = labels[i]
        if label is None:
            label = os.path.splitext(fname)[0]
            label = label.split('_fields')[0]

        plt.plot(t, vrms, label=label, linewidth=2, color=colors[i])

    plt.xlabel(r'time (s)')
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    plt.ylabel(r'rms volume error')
    plt.grid(linestyle='dashed', zorder=-1)
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('rms_vol_err.' + fmt, bbox_inches='tight')
    plt.close()


def plot_max_volume_error(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the gridded absolute volume error (normalised with
    cell volume).
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    plt.figure()

    h5reader = H5Reader()
    for i, fname in enumerate(fnames):
        h5reader.open(fname)

        if h5reader.is_parcel_file:
            raise IOError('Not a field output file.')

        vmax = h5reader.get_diagnostic('max absolute normalised volume error')
        t = h5reader.get_diagnostic('t')
        h5reader.close()

        label = labels[i]
        if label is None:
            label = os.path.splitext(fname)[0]
            label = label.split('_fields')[0]

        plt.plot(t, vmax, label=label, linewidth=2, color=colors[i])

    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    plt.xlabel(r'time (s)')
    plt.ylabel(r'max normalised volume error')
    plt.grid(linestyle='dashed', zorder=-1)
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))
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

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    lam_mean = np.zeros(nsteps)
    lam_std = np.zeros(nsteps)
    t = np.zeros(nsteps)

    for step in range(nsteps):
        lam = h5reader.get_aspect_ratio(step)

        lam_mean[step] = lam.mean()
        lam_std[step] = lam.std()

        t[step] = h5reader.get_step_attribute(step, 't')

    lmax = h5reader.get_parcel_option('lambda')[0]

    h5reader.close()

    plt.figure()
    plt.plot(t, lam_mean, color='blue', label=r'mean')
    plt.fill_between(t, lam_mean - lam_std, lam_mean + lam_std,
                     alpha=0.5, label=r'std. dev.')

    plt.axhline(lmax, linestyle='dashed', color='black',
                label=r'$\lambda\le\lambda_{\max} = ' + str(lmax) + '$')

    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'time (s)')
    plt.ylabel(r'aspect ratio $\lambda$')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))

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

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    nparcels = np.zeros(nsteps)
    t = np.zeros(nsteps)

    for step in range(nsteps):
        nparcels[step] = h5reader.get_num_parcels(step)
        t[step] = h5reader.get_step_attribute(step, 't')

    h5reader.close()

    plt.figure()
    plt.plot(t, nparcels, color='blue')
    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'time (s)')
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

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    vol_mean = np.zeros(nsteps)
    vol_std = np.zeros(nsteps)
    t = np.zeros(nsteps)

    extent = h5reader.get_box_extent()
    ncells = h5reader.get_box_ncells()
    vcell = np.prod(extent / ncells)

    for step in range(nsteps):
        vol = h5reader.get_dataset(step, 'volume')

        vol_mean[step] = vol.mean() / vcell
        vol_std[step] = vol.std() / vcell

        t[step] = h5reader.get_step_attribute(step, 't')

    h5reader.close()

    plt.figure()
    plt.plot(t, vol_mean, color='blue', label=r'mean')
    plt.fill_between(t, vol_mean - vol_std, vol_mean + vol_std,
                     alpha=0.5, label=r'std. dev.')

    plt.axhline(1.0, linestyle='dashed', color='black',
                label=r'cell volume $V_{0}$')

    plt.grid(linestyle='dashed', zorder=-1)

    plt.xlabel(r'time (s)')
    plt.ylabel(r'parcel volume / $V_{0}$')

    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15))

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_parcel_volume_profile.' + fmt, bbox_inches='tight')
    plt.close()


def plot_center_of_mass(fname, show=False, fmt="png"):
    prefix, ext = os.path.splitext(fname)

    if ext == '.hdf5':
        print('Extract data, evaluate quantities, write CSV file and plot.')

        h5reader = H5Reader()
        h5reader.open(fname)

        if not h5reader.is_parcel_file:
            raise IOError('Not a parcel output file.')

        nsteps = h5reader.get_num_steps()

        bi_vi = np.zeros(nsteps-1)
        bi_vi_xi = np.zeros(nsteps-1)
        bi_vi_zi = np.zeros(nsteps-1)

        bi_vi_x2i = np.zeros(nsteps-1)
        bi_vi_z2i = np.zeros(nsteps-1)

        wi_vi = np.zeros(nsteps-1)
        wi_vi_xi = np.zeros(nsteps-1)
        wi_vi_zi = np.zeros(nsteps-1)

        wi_vi_x2i = np.zeros(nsteps-1)
        wi_vi_z2i = np.zeros(nsteps-1)

        xi_zi = np.zeros(nsteps-1)

        t = np.zeros(nsteps-1)

        # skip zero step since vorticity is not given
        for j in range(1, nsteps):

            pos = h5reader.get_dataset(j, 'position')
            vol = h5reader.get_dataset(j, 'volume')
            vor = h5reader.get_dataset(j, 'vorticity')
            buo = h5reader.get_dataset(j, 'buoyancy')
            t[j-1] = h5reader.get_step_attribute(j, 't')

            # we only want parcels with x > 0
            ind = (pos[0, :] > 0.0)

            pos = pos[:, ind]
            vol = vol[ind]
            vor = vor[ind]
            buo = buo[ind]

            bv = buo * vol
            vv = vor * vol

            # first moment
            bi_vi[j-1] = bv.mean()
            bi_vi_xi[j-1] = (bv * pos[0, :]).mean()
            bi_vi_zi[j-1] = (bv * pos[1, :]).mean()

            # second moment
            bi_vi_x2i[j-1] = (bv * pos[0, :] ** 2).mean()
            bi_vi_z2i[j-1] = (bv * pos[1, :] ** 2).mean()

            # first moment
            wi_vi[j-1] = vv.mean()
            wi_vi_xi[j-1] = (vv * pos[0, :]).mean()
            wi_vi_zi[j-1] = (vv * pos[1, :]).mean()

            # second moment
            wi_vi_x2i[j-1] = (vv * pos[0, :] ** 2).mean()
            wi_vi_z2i[j-1] = (vv * pos[1, :] ** 2).mean()

            bi_vi_xi_zi = (bv * pos[0, :] * pos[1, :]).mean()
            wi_vi_xi_zi = (vv * pos[0, :] * pos[1, :]).mean()

        h5reader.close()

        xb_bar = bi_vi_xi / bi_vi
        zb_bar = bi_vi_zi / bi_vi
        xw_bar = wi_vi_xi / wi_vi
        zw_bar = wi_vi_zi / wi_vi

        data = {
            't':        t,
            'xb_bar':   xb_bar,
            'zb_bar':   zb_bar,
            'x2b_bar':  bi_vi_x2i / bi_vi - xb_bar ** 2,
            'z2b_bar':  bi_vi_z2i / bi_vi - zb_bar ** 2,
            'xzb_bar':  bi_vi_xi_zi / bi_vi - xb_bar * zb_bar,
            'xw_bar':   xw_bar,
            'zw_bar':   zw_bar,
            'x2w_bar':  wi_vi_x2i / wi_vi - xw_bar ** 2,
            'z2w_bar':  wi_vi_z2i / wi_vi - zw_bar ** 2,
            'xzw_bar':  wi_vi_xi_zi / wi_vi - xw_bar * zw_bar
        }

        df = pd.DataFrame(data=data)

        df.to_csv(prefix + '.csv', index=False)
    elif ext == '.csv':
        print('Read CSV file and plot.')
        df = pd.read_csv(prefix + '.csv')
    else:
        raise IOError('Wrong file format. Requires parcel hdf5 or CSV file.')


    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(9, 6), dpi=200)
    axes[0].plot(df['t'], df['xb_bar'], label=r'$\langle x\rangle_b$')

    std = np.sqrt(df['x2b_bar'])
    axes[0].fill_between(df['t'], df['xb_bar']-std, df['xb_bar']+std, alpha=0.5,
                         label=r'$\pm\sqrt{\langle x^2\rangle_b}$')

    axes[0].plot(df['t'], df['zb_bar'], label=r'$\langle z\rangle_b$')

    std = np.sqrt(df['z2b_bar'])
    axes[0].fill_between(df['t'], df['zb_bar']-std, df['zb_bar']+std, alpha=0.5,
                         label=r'$\pm\sqrt{\langle z^2\rangle_b}$')


    axes[1].plot(df['t'], df['xw_bar'], label=r'$\langle x\rangle_\zeta$')

    std = np.sqrt(df['x2w_bar'])
    axes[1].fill_between(df['t'], df['xw_bar']-std, df['xw_bar']+std, alpha=0.5,
                         label=r'$\pm\sqrt{\langle x^2\rangle_\zeta}$')

    axes[1].plot(df['t'], df['zw_bar'], label=r'$\langle z\rangle_\zeta$')

    std = np.sqrt(df['z2w_bar'])
    axes[1].fill_between(df['t'], df['zw_bar']-std, df['zw_bar']+std, alpha=0.5,
                         label=r'$\pm\sqrt{\langle z^2\langle_\zeta}$')

    axes[0].grid(which='both', linestyle='dashed')
    axes[0].legend(loc='upper right', ncol=1, bbox_to_anchor=(1.25, 1.0))

    axes[1].grid(which='both', linestyle='dashed')
    axes[1].legend(loc='upper right', ncol=1, bbox_to_anchor=(1.25, 1.0))
    axes[1].set_xlabel(r'time (s)')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig(prefix + '_center_of_mass.' + fmt,
                    bbox_inches='tight')
    plt.close()



def plot_cumulative(fnames, step=0, dset='volume', show=False, fmt="png", **kwargs):
    """
    Plot the mean and standard deviation of the parcel volume
    normalised with the cell volume.
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    for i, fname in enumerate(fnames):
        h5reader = H5Reader()
        h5reader.open(fname)

        if not h5reader.is_parcel_file:
            raise IOError('Not a parcel output file.')

        nsteps = h5reader.get_num_steps()

        if step < 0 or step > nsteps-1:
            raise ValueError('Number of steps exceeded.')


        data = {dset: h5reader.get_dataset(step, dset)}

        h5reader.close()

        label = labels[i]
        if label is None:
            label = os.path.splitext(fname)[0]
            label = label.split('_parcels')[0]

        sns.ecdfplot(data=data, x=dset, stat='proportion',
                     color=colors[i], label=label)

    plt.ylabel('proportion')
    plt.grid(which='both', linestyle='dashed')

    if n > 1:
        plt.legend(loc='upper center', ncol=min(4, n),
                   bbox_to_anchor=(0.5, 1.15))

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        plt.savefig(prefix + 'parcel_cumulative_' + dset + '.' + fmt,
                    bbox_inches='tight')
    plt.close()


def plot_time_pie_chart(fname, show=False, fmt="png"):
    df = pd.read_csv(fname)

    # remove epic
    epic_time = df['total time'][0]
    df.drop([0], inplace=True)

    ind = df['percentage'] < 1.0
    # 25 June 2021
    # https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
    df.drop(df[ind].index, inplace=True)

    others = epic_time - df['total time'].sum()

    df2 = pd.DataFrame([['others', 0, others, others / epic_time * 100.0]],
                       columns=df.columns)

    df = df.append(df2)

    df.sort_values(by=['total time'], inplace=True, ascending=False)

    n = len(df['function name'])

    plt.figure(figsize=(7, 5))
    ax = plt.gca()

    data = list(df['total time'])

    # 25 June 2021
    # https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_and_donut_labels.html#sphx-glr-gallery-pie-and-polar-charts-pie-and-donut-labels-py
    wedges, texts, autotexts = ax.pie(data, autopct=lambda pct: get_autopct(pct, data), pctdistance=0.8,
                                      wedgeprops=dict(width=0.5), startangle=30,
                                      textprops=dict(color="w"))

    bbox_props = dict(boxstyle="round,pad=0.3",
                      facecolor='wheat', alpha=0.5)
    kw = dict(arrowprops=dict(arrowstyle="-"),
             bbox=bbox_props, zorder=0, va="center")

    labels = list(df['function name'])

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

    plt.setp(autotexts, size=8, weight="bold")
    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fname)[0]
        plt.savefig(prefix + '_timing_pie_chart.' + fmt,
                    bbox_inches='tight')
    plt.close()

def plot_time_bar(fnames, show=False, fmt="png", **kwargs):

    labels = kwargs.pop('labels', len(fnames) * [None])

    if len(labels) < len(fnames):
        raise ValueError('Not enough labels provided.')

    df = pd.read_csv(fnames[0])

    # remove epic timer
    epic_time = df['total time'][0]
    df.drop([0], inplace=True)

    ind = df['percentage'] < 1.0
    df.drop(df[ind].index, inplace=True)

    # remove unncessary columns
    df.drop(['#calls', 'percentage'], axis=1, inplace=True)

    others = epic_time - df['total time'].sum()

    df2 = pd.DataFrame([['others', others]], columns=df.columns)
    df = df.append(df2)
    df.sort_values(by=['total time'], inplace=True, ascending=False)

    # rename column
    label = labels[0]
    if label is None:
        label = os.path.splitext(fnames[0])[0]

    # 27 June 2021
    # https://stackoverflow.com/questions/11346283/renaming-columns-in-pandas
    df.columns = ['function name', label]

    names = df['function name'].copy()


    for n in range(1, len(fnames)):
        df2 = pd.read_csv(fnames[n])

        epic_time = df2['total time'][0]
        df2.drop([0], inplace=True)

        df2.drop(df2[ind].index, inplace=True)

        # remove unncessary columns
        df2.drop(['#calls', 'percentage'], axis=1, inplace=True)

        others = epic_time - df2['total time'].sum()

        df3 = pd.DataFrame([['others', others]], columns=df2.columns)
        df2 = df2.append(df3)

        # https://stackoverflow.com/questions/26707171/sort-pandas-dataframe-based-on-list
        # 27 June 2021
        df2['tmp'] = pd.Categorical(df2['function name'], categories=names, ordered=True)
        df2.sort_values('tmp', inplace=True)

        # rename column
        label = labels[n]
        if label is None:
            label = os.path.splitext(fnames[n])[0]

        # 27 June 2021
        # https://stackoverflow.com/questions/27965295/dropping-rows-from-dataframe-based-on-a-not-in-condition
        df2 = df2[df2['function name'].isin(names)]

        df2.drop(['function name', 'tmp'], axis=1, inplace=True)
        df2.columns = [label]

        # append to df
        df = pd.concat([df, df2], axis=1)


    df.plot.bar(x='function name')

    # 26 June 2021
    # https://stackoverflow.com/questions/14852821/aligning-rotated-xticklabels-with-their-respective-xticks
    plt.xticks(rotation=30, ha='right')

    plt.ylabel(r'wall time (s)')
    plt.xlabel('')
    plt.grid(which='both', linestyle='dashed')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('timing_bar_plot.' + fmt, bbox_inches='tight')
    plt.close()


#def plot_field(fname, show=False, fmt="png", ax=None):

    #h5reader = H5Reader()
    #h5reader.open(fname)

    #extent = h5reader.get_box_extent()
    #origin = h5reader.get_box_origin()
    #ncells = h5reader.get_box_ncells()

    #xx = np.linspace(origin[0], origin[0] + extent[0], ncells[0], endpoint=False)
    #yy = np.linspace(origin[1], origin[1] + extent[1], ncells[1], endpoint=False)
    #yy, xx = np.boxgrid(yy, xx)

#sc = axes.pcolorbox(xx, yy, npc, shading='nearest')
#plt.colorbar(sc, ax=axes)
#plt.xlabel('x')
#plt.ylabel('y')
#plt.tight_layout()
#plt.show()
#plt.close()
