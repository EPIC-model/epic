from tools.h5_reader import H5Reader
from tools.plot_beautify import *
from tools.plot_style import *
from tools.units import *
import matplotlib.colors as cls
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns

def _plot_parcels(ax, h5reader, step, coloring, vmin, vmax, draw_cbar=True, **kwargs):

    # 19 Feb 2021
    # https://stackoverflow.com/questions/43009724/how-can-i-convert-numbers-to-a-color-scale-in-matplotlib
    norm = cls.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis_r

    origin = h5reader.get_box_origin()
    extent = h5reader.get_box_extent()

    # instantiating the figure object
    fkwargs = {k: v for k, v in kwargs.items() if v is not None}
    left = fkwargs.get('xmin', origin[0])
    right = fkwargs.get('xmax', origin[0] + extent[0])
    bottom = fkwargs.get('ymin', origin[1])
    top = fkwargs.get('ymax', origin[1] + extent[1])

    if coloring == 'aspect-ratio':
        data = h5reader.get_aspect_ratio(step=step)
    elif coloring == 'vol-distr':
        data = h5reader.get_dataset(step=step, name='volume')
        # 5 August 2021
        # https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        # https://stackoverflow.com/questions/40601997/setting-discrete-colormap-corresponding-to-specific-data-range-in-matplotlib
        cmap = plt.cm.get_cmap('bwr', 2)
        bounds = [0, vmin, vmax]
        norm = cls.BoundaryNorm(bounds, cmap.N)
    else:
        data = h5reader.get_dataset(step=step, name=coloring)

    ells = h5reader.get_ellipses(step=step)
    for j, e in enumerate(ells):
        ax.add_artist(e)
        #e.set_clip_box(ax.bbox)
        e.set_alpha(0.75)
        e.set_facecolor(cmap(norm(data[j])))

    ax.set_xlim([left, right])
    ax.set_ylim([bottom, top])

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

        if coloring == 'aspect-ratio':
            cbar.set_label(r'$1 \leq \lambda \leq \lambda_{max}$')
        elif coloring == 'vol-distr':
            # 5 August 2021
            # https://matplotlib.org/stable/gallery/ticks_and_spines/colorbar_tick_labelling_demo.html
            cbar.ax.set_yticklabels([r'0', r'$V_{min}$', r'$V_{max}$'])
        else:
            cbar.set_label(coloring)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')


def plot_parcels(fname, step, show=False, fmt="png",
                  coloring='aspect-ratio', **kwargs):
    h5reader = H5Reader()

    h5reader.open(fname)

    if not h5reader.is_parcel_file:
        raise IOError('Not a parcel output file.')

    nsteps = h5reader.get_num_steps()

    if step > nsteps - 1:
        raise ValueError('Step number exceeds limit of ' + str(nsteps-1) + '.')

    if step < 0:
        raise ValueError('Step number cannot be negative.')

    if coloring == 'aspect-ratio':
        vmin = 1.0
        vmax = h5reader.get_parcel_option('lambda')
    elif coloring == 'vol-distr':
        extent = h5reader.get_box_extent()
        ncells = h5reader.get_box_ncells()
        vcell = np.prod(extent / ncells)
        vmin = vcell / h5reader.get_parcel_option('vmin_fraction')
        vmax = vcell / h5reader.get_parcel_option('vmax_fraction')
    else:
        vmin, vmax = h5reader.get_dataset_min_max(coloring)

    plt.figure(num=step)

    _plot_parcels(plt.gca(), h5reader, step, coloring, vmin, vmax, **kwargs)

    if show:
        plt.tight_layout()
        plt.show()
    else:
        plt.savefig('parcels_'  + coloring + '_step_' + str(step).zfill(len(str(nsteps))) + '.' + fmt,
                    bbox_inches='tight')
    plt.close()
    h5reader.close()


def plot_volume_symmetry_error(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the symmetry error of the gridded volume.
    (The gridded symmetry volume is only written in debug mode.)
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    plt.figure()
    for i, fname in enumerate(fnames):
        h5reader = H5Reader()
        h5reader.open(fname)

        if h5reader.is_parcel_file:
            raise IOError('Not a field output file.')

        try:
            h5reader.get_step_attribute(0, 'max symmetry volume error')
        except:
            raise IOError('This plot is only available in debug mode.')

        nsteps = h5reader.get_num_steps()

        vmax = np.zeros(nsteps)
        t = np.zeros(nsteps)
        for step in range(nsteps):
            vmax[step] = h5reader.get_step_attribute(step, 'max symmetry volume error')
            t[step] = h5reader.get_step_attribute(step, 't')

        h5reader.close()

        plt.fill_between(t, 0, vmax, color=colors[i], label=labels[i],
                        edgecolor=colors[i], linewidth=0.75)

    plt.grid(which='both', linestyle='dashed')
    plt.xlabel(get_label('time', units['time']))
    plt.ylabel(r'volume symmetry error')
    plt.yscale('log')
    #plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))


    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        plt.savefig(prefix + 'vol_sym_err.' + fmt,
                    bbox_inches='tight')
    plt.close()


def plot_rms_volume_error(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the gridded rms volume error.
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])
    yscale = kwargs.pop('yscale', 'linear')
    ylim = kwargs.pop('ylim', (None, None))

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

        plt.plot(t, vrms, label=labels[i], linewidth=2, color=colors[i])

    plt.xlabel(get_label('time', units['time']))
    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    plt.ylabel(r'rms area error')
    plt.grid(which='both', linestyle='dashed', zorder=-1)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])
    plt.yscale(yscale)
    plt.ylim(ylim)
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

        plt.plot(t, vmax, label=labels[i], linewidth=2, color=colors[i])

    plt.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    plt.xlabel(get_label('time', units['time']))
    plt.ylabel(r'max normalised volume error')
    plt.grid(linestyle='dashed', zorder=-1)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('max_normalised_vol_err.' + fmt, bbox_inches='tight')
    plt.close()


def plot_parcel_profile(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the mean and standard deviation of the parcel aspect ratio.
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])
    dset = kwargs.pop('dset', 'aspect-ratio')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    h5reader = H5Reader()

    plt.figure()

    lmax = 0

    for i, fname in enumerate(fnames):

        h5reader.open(fname)

        if not h5reader.is_parcel_file:
            raise IOError('Not a parcel output file.')

        nsteps = h5reader.get_num_steps()

        data_mean = np.zeros(nsteps)
        data_std = np.zeros(nsteps)
        t = np.zeros(nsteps)

        for step in range(nsteps):

            data = None

            if dset == 'aspect-ratio':
                data = h5reader.get_aspect_ratio(step)
            else:
                data = h5reader.get_dataset(step, dset)


            if dset == 'volume':
                extent = h5reader.get_box_extent()
                ncells = h5reader.get_box_ncells()
                vcell = np.prod(extent / ncells)
                data /= vcell


            data_mean[step] = data.mean()
            data_std[step] = data.std()

            t[step] = h5reader.get_step_attribute(step, 't')

        if dset == 'aspect-ratio':
            lmax = max(lmax, h5reader.get_parcel_option('lambda'))

        h5reader.close()

        plt.plot(t, data_mean, label=labels[i], color=colors[i])
        plt.fill_between(t, data_mean - data_std, data_mean + data_std,
                         alpha=0.5, color=colors[i])


    plt.xlabel(get_label('time', units['time']))
    plt.grid(linestyle='dashed', zorder=-1)

    if dset == 'aspect-ratio':
        plt.ylabel(r'aspect ratio $\lambda$')
        plt.text(t[10], 0.92 * lmax, r'$\lambda\le\lambda_{max} = ' + str(lmax) + '$')
        plt.axhline(lmax, linestyle='dashed', color='black')
    elif dset == 'volume':
        plt.ylabel(r'parcel volume / $V_{g}$')
        #plt.axhline(1.0, linestyle='dashed', color='black',
                #label=r'cell volume $V_{g}$')
    else:
        plt.ylabel(r'parcel ' + dset)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        plt.savefig(prefix + 'parcel_' + dset + '_profile.' + fmt,
                    bbox_inches='tight')
    plt.close()


def plot_parcel_stats_profile(fnames, show=False, fmt="png", **kwargs):
    """
    Plot parcel statistics
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])
    dset = kwargs.pop('dset', 'aspect-ratio')

    if dset == 'aspect-ratio':
        dset = 'aspect ratio'

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    h5reader = H5Reader()

    plt.figure()

    lmax = 0

    for i, fname in enumerate(fnames):

        h5reader.open(fname)

        if not h5reader.is_parcel_stats_file:
            raise IOError('Not a parcel diagnostic output file.')

        nsteps = h5reader.get_num_steps()

        data_mean = np.zeros(nsteps)
        data_std = np.zeros(nsteps)
        t = np.zeros(nsteps)

        for step in range(nsteps):
            t[step] = h5reader.get_step_attribute(step, 't')
            data_mean[step] = h5reader.get_step_attribute(step, 'avg ' + dset)
            data_std[step] = h5reader.get_step_attribute(step, 'std ' + dset)

        if dset == 'aspect ratio':
            lmax = max(lmax, h5reader.get_parcel_option('lambda'))

        h5reader.close()

        plt.plot(t, data_mean, label=labels[i], color=colors[i])
        plt.fill_between(t, data_mean - data_std, data_mean + data_std,
                         alpha=0.5, color=colors[i])


    plt.xlabel(get_label('time', units['time']))
    plt.grid(linestyle='dashed', zorder=-1)

    if dset == 'aspect-ratio':
        plt.ylabel(r'aspect ratio $\lambda$')
        plt.text(t[10], 0.95 * lmax, r'$\lambda\le\lambda_{max} = ' + str(lmax) + '$')
        plt.axhline(lmax, linestyle='dashed', color='black')
    elif dset == 'volume':
        plt.ylabel(r'parcel volume / $V_{g}$')
        #plt.axhline(1.0, linestyle='dashed', color='black',
                #label=r'cell volume $V_{g}$')
    else:
        plt.ylabel(r'parcel ' + dset)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        dset = dset.replace(' ', '_')
        plt.savefig(prefix + 'parcel_' + dset + '_profile.' + fmt,
                    bbox_inches='tight')
    plt.close()


def plot_parcel_number(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the number of parcels in simulation.
    """
    labels = kwargs.pop('labels', None)

    if labels is None:
        labels = [None] * len(fnames)

    plt.figure()
    for i, fname in enumerate(fnames):

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

        plt.plot(t, nparcels, label=labels[i])

    plt.grid(linestyle='dashed', zorder=-1)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'],
                   ncol=min(len(labels), legend_dict['ncol']),
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.xlabel(get_label('time', units['time']))
    plt.ylabel(r'parcel count')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('parcel_number_profile.' + fmt, bbox_inches='tight')
    plt.close()


def plot_small_parcel_number(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the number of small parcels in simulation.
    """
    labels = kwargs.pop('labels', None)

    if labels is None:
        labels = [None] * len(fnames)

    plt.figure()
    for i, fname in enumerate(fnames):

        h5reader = H5Reader()
        h5reader.open(fname)

        if not h5reader.is_parcel_stats_file:
            raise IOError('Not a parcel diagnostic output file.')

        nsteps = h5reader.get_num_steps()

        nparcels = np.zeros(nsteps)
        nsmall = np.zeros(nsteps)
        t = np.zeros(nsteps)

        for step in range(nsteps):
            nparcels[step] = h5reader.get_step_attribute(step, 'num parcel')
            nsmall[step] = h5reader.get_step_attribute(step, 'num small parcels')
            t[step] = h5reader.get_step_attribute(step, 't')

        h5reader.close()

        plt.plot(t, nsmall / nparcels * 100.0, label=labels[i])

    plt.grid(linestyle='dashed', zorder=-1)

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'],
                   ncol=min(len(labels), legend_dict['ncol']),
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.xlabel(get_label('time', units['time']))
    plt.ylabel(r'percentage of small parcels (\%)')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('parcel_small_number_profile.' + fmt, bbox_inches='tight')
    plt.close()

def plot_center_of_mass(fnames, show=False, fmt="png", dset='buoyancy', **kwargs):

    tag = None
    if dset == 'buoyancy':
        tag = 'b'
        tag_name = 'b'
    elif dset == 'vorticity':
        tag = 'w'
        tag_name = '\zeta'
    else:
        raise ValueError("Dataset must be 'buoyancy' or 'vorticity'.")

    labels = kwargs.pop('labels', None)
    variance = kwargs.pop('variance', False)

    n = len(fnames)
    if labels is None:
        labels = [None] * n

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    h5reader = H5Reader()

    fig1 = plt.figure(num=1)
    ax1 = fig1.gca()

    fig2 = plt.figure(num=2)
    ax2 = fig2.gca()

    for i, fname in enumerate(fnames):

        h5reader.open(fname)

        if not h5reader.is_parcel_stats_file:
            raise IOError('Not a parcel diagnostic output file.')

        nsteps = h5reader.get_num_steps()

        t = np.zeros(nsteps)

        xb_bar = np.zeros(nsteps)
        zb_bar = np.zeros(nsteps)
        x2b_bar = np.zeros(nsteps)
        z2b_bar = np.zeros(nsteps)
        xzb_bar = np.zeros(nsteps)
        xw_bar = np.zeros(nsteps)
        zw_bar = np.zeros(nsteps)
        x2w_bar = np.zeros(nsteps)
        z2w_bar = np.zeros(nsteps)
        xzw_bar = np.zeros(nsteps)


        ## skip zero step since vorticity is not given
        for j in range(0, nsteps):

            t[j] = h5reader.get_step_attribute(j, 't')
            xb_bar[j] = h5reader.get_step_attribute(j, 'xb_bar')
            zb_bar[j] = h5reader.get_step_attribute(j, 'zb_bar')
            x2b_bar[j] = h5reader.get_step_attribute(j, 'x2b_bar')
            z2b_bar[j] = h5reader.get_step_attribute(j, 'z2b_bar')
            xzb_bar[j] = h5reader.get_step_attribute(j, 'xzb_bar')
            xw_bar[j] = h5reader.get_step_attribute(j, 'xv_bar')
            zw_bar[j] = h5reader.get_step_attribute(j, 'zv_bar')
            x2w_bar[j] = h5reader.get_step_attribute(j, 'x2v_bar')
            z2w_bar[j] = h5reader.get_step_attribute(j, 'z2v_bar')
            xzw_bar[j] = h5reader.get_step_attribute(j, 'xzv_bar')


        h5reader.close()

        data = {
            't':        t,
            'xb_bar':   xb_bar,
            'zb_bar':   zb_bar,
            'x2b_bar':  x2b_bar,
            'z2b_bar':  z2b_bar,
            'xzb_bar':  xzb_bar,
            'xw_bar':   xw_bar,
            'zw_bar':   zw_bar,
            'x2w_bar':  x2w_bar,
            'z2w_bar':  z2w_bar,
            'xzw_bar':  xzw_bar,
        }

        df = pd.DataFrame(data=data)

        ax1.plot(df['t'], df['x' + tag + '_bar'], label=labels[i], color=colors[i])

        if variance:
            std = np.sqrt(df['x2' + tag + '_bar'])
            ax1.fill_between(df['t'], df['x' + tag + '_bar']-std, df['x' + tag + '_bar']+std, alpha=0.5,
                             color=colors[i], edgecolor='None')

        ax2.plot(df['t'], df['z' + tag + '_bar'], label=labels[i], color=colors[i])

        if variance:
            std = np.sqrt(df['z2' + tag + '_bar'])
            ax2.fill_between(df['t'], df['z' + tag + '_bar']-std, df['z' + tag + '_bar']+std, alpha=0.5,
                             color=colors[i], edgecolor='None')

    ax1.grid(which='both', linestyle='dashed')
    ax2.grid(which='both', linestyle='dashed')

    if not labels[0] is None:
        ax1.legend(loc=legend_dict['loc'],
                   ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])
        ax2.legend(loc=legend_dict['loc'],
                   ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    ax1.set_xlabel(r'time (s)')
    ax2.set_xlabel(r'time (s)')
    if variance:
        ax1.set_ylabel(r'$\langle x\rangle_' + tag_name + \
            r'\pm\sqrt{\langle x^2\rangle_' + tag_name + '}$')
        ax2.set_ylabel(r'$\langle y\rangle_' + tag_name + \
            r'\pm\sqrt{\langle y^2\rangle_' + tag_name + '}$')
    else:
        ax1.set_ylabel(r'$\langle x\rangle_' + tag_name + r'$')
        ax2.set_ylabel(r'$\langle y\rangle_' + tag_name + r'$')
    fig1.tight_layout()
    fig2.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        fig1.savefig(prefix + 'x_center_of_mass_' + dset + '.' + fmt,
                     bbox_inches='tight')
        fig2.savefig(prefix + 'z_center_of_mass_' + dset + '.' + fmt,
                     bbox_inches='tight')
    plt.close('all')


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

        sns.ecdfplot(data=data, x=dset, stat='proportion',
                     color=colors[i], label=labels[i])

    plt.ylabel('proportion')
    plt.grid(which='both', linestyle='dashed')

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

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


def plot_parcels_per_cell(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the mean and standard deviation of the number of
    parcels per cell
    """
    n = len(fnames)

    labels = kwargs.pop('labels', n * [None])
    add_minmax = kwargs.pop('add_minmax', True)

    if len(labels) < n:
        raise ValueError('Not enough labels provided.')

    colors =  plt.cm.tab10(np.arange(n).astype(int))

    for i, fname in enumerate(fnames):
        h5reader = H5Reader()
        h5reader.open(fname)

        if not h5reader.is_field_stats_file:
            raise IOError('Not a field diagnostic output file.')

        nsteps = h5reader.get_num_steps()

        n_avg = np.zeros(nsteps)

        if add_minmax:
            n_min = np.zeros(nsteps)
            n_max = np.zeros(nsteps)

        t = np.zeros(nsteps)

        for step in range(nsteps):
            n_avg[step] = h5reader.get_step_attribute(step, 'average num parcels per cell')
            t[step] = h5reader.get_step_attribute(step, 't')

            if add_minmax:
                n_min[step] = h5reader.get_step_attribute(step, 'min num parcels per cell')
                n_max[step] = h5reader.get_step_attribute(step, 'max num parcels per cell')

        h5reader.close()

        plt.plot(t, n_avg, color=colors[i], label=labels[i])

        if add_minmax:
            plt.fill_between(t, n_min, n_max, color=colors[i], alpha=0.5, edgecolor=None)

    plt.xlabel(get_label('time', units['time']))
    plt.ylabel(r'number of parcels/cell')
    plt.grid(which='both', linestyle='dashed')

    if not labels[0] is None:
        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

    plt.tight_layout()

    if show:
        plt.show()
    else:
        prefix = os.path.splitext(fnames[0])[0] + '_'
        if n > 1:
            prefix = ''
        plt.savefig(prefix + 'number_of_parcels_per_cell.' + fmt,
                    bbox_inches='tight')
    plt.close()


def plot_energy(fnames, show=False, fmt="png", **kwargs):
    """
    Plot the kinetic, potential and total energy.
    """
    n = len(fnames)

    for i, fname in enumerate(fnames):

        plt.figure(num=i)

        h5reader = H5Reader()
        h5reader.open(fname)

        if not h5reader.is_parcel_stats_file:
            raise IOError('Not a parcel diagnostic output file.')

        nsteps = h5reader.get_num_steps()

        pe = np.zeros(nsteps)
        ke = np.zeros(nsteps)
        te = np.zeros(nsteps)
        t = np.zeros(nsteps)

        for step in range(nsteps):
            pe[step] = h5reader.get_step_attribute(step, 'potential energy')
            ke[step] = h5reader.get_step_attribute(step, 'kinetic energy')
            te[step] = h5reader.get_step_attribute(step, 'total energy')
            t[step] = h5reader.get_step_attribute(step, 't')

        h5reader.close()

        plt.plot(t, pe, label=r'$\mathcal{P}$')
        plt.plot(t, ke, label=r'$\mathcal{K}$')
        plt.plot(t, te, label=r'$\mathcal{P}+\mathcal{K}$')
        plt.xlabel(get_label('time', units['time']))
        plt.ylabel(r'energy')
        plt.grid(which='both', linestyle='dashed')

        plt.legend(loc=legend_dict['loc'], ncol=legend_dict['ncol'],
                   bbox_to_anchor=legend_dict['bbox_to_anchor'])

        plt.tight_layout()

        if show:
            plt.show()
        else:
            prefix = os.path.splitext(fname)[0]
            prefix = prefix.split('parcel_stats')[0]
            plt.savefig(prefix + 'energy.' + fmt,
                        bbox_inches='tight')
        plt.close()


def plot_time_pie_chart(fname, show=False, fmt="png", **kwargs):
    df = pd.read_csv(fname)

    # remove epic
    epic_time = df['total time'][0]
    df.drop([0], inplace=True)

    skip = kwargs.pop('skip', 1.0)

    ind = df['percentage'] < skip
    # 25 June 2021
    # https://stackoverflow.com/questions/13851535/how-to-delete-rows-from-a-pandas-dataframe-based-on-a-conditional-expression
    df.drop(df[ind].index, inplace=True)

    others = epic_time - df['total time'].sum()

    df2 = pd.DataFrame([['others', 0, others, others / epic_time * 100.0]],
                       columns=df.columns)

    df = df.append(df2)

    df.sort_values(by=['total time'], inplace=True, ascending=False)

    n = len(df['function name'])

    plt.figure()
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
    skip = kwargs.pop('skip', 1.0)

    if len(labels) < len(fnames):
        raise ValueError('Not enough labels provided.')

    df = pd.read_csv(fnames[0])

    # remove epic timer
    epic_time = df['total time'][0]
    df.drop([0], inplace=True)

    ind = df['percentage'] < skip
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


def plot_time_speedup(fnames, nthreads, show=False, fmt="png"):

    if len(fnames) < 2:
        raise ValueError('Not enough files provided.')

    if not len(fnames) == len(nthreads):
        raise RuntimeError('List of files not equal to thread list.')

    df = pd.read_csv(fnames[0])

    epic_time = df['total time'][0]

    ind = df['percentage'] < 1.0
    df.drop(df[ind].index, inplace=True)

    # remove unncessary columns
    df.drop(['#calls', 'percentage'], axis=1, inplace=True)

    # multiply by 2 since EPIC timer is not removed from df
    others = 2.0 * epic_time - df['total time'].sum()

    df2 = pd.DataFrame([['others', others]], columns=df.columns)
    df = df.append(df2)
    df.sort_values(by=['total time'], inplace=True, ascending=False)

    # rename column
    # 27 June 2021
    # https://stackoverflow.com/questions/11346283/renaming-columns-in-pandas
    df.columns = ['function name', '0']

    names = df['function name'].copy()


    for n in range(1, len(fnames)):
        df2 = pd.read_csv(fnames[n])

        epic_time = df2['total time'][0]

        df2.drop(df2[ind].index, inplace=True)

        # remove unncessary columns
        df2.drop(['#calls', 'percentage'], axis=1, inplace=True)

        # multiply by 2 since EPIC timer is not removed from df
        others = 2.0 * epic_time - df2['total time'].sum()

        df3 = pd.DataFrame([['others', others]], columns=df2.columns)
        df2 = df2.append(df3)

        # https://stackoverflow.com/questions/26707171/sort-pandas-dataframe-based-on-list
        # 27 June 2021
        df2['tmp'] = pd.Categorical(df2['function name'], categories=names, ordered=True)
        df2.sort_values('tmp', inplace=True)

        # rename column
        label = str(n)

        # 27 June 2021
        # https://stackoverflow.com/questions/27965295/dropping-rows-from-dataframe-based-on-a-not-in-condition
        df2 = df2[df2['function name'].isin(names)]

        df2.drop(['function name', 'tmp'], axis=1, inplace=True)
        df2.columns = [label]

        # append to df
        df = pd.concat([df, df2], axis=1)


    df.reset_index(drop=True, inplace=True)
    names = list(df['function name'])
    df.reindex(names)
    df.drop('function name', axis=1, inplace=True)

    for i, name in enumerate(names):
        # get i-th row
        data = df.loc[i, :]
        plt.semilogy(nthreads, data, label=name)

    # ideal timing
    data = df.loc[0, :][0]
    data = np.asarray(len(nthreads) * [data]) / np.asarray(nthreads)
    plt.semilogy(nthreads, data, linestyle='dashed', color='black', label='perfect scaling')

    plt.ylabel(r'wall time (s)')
    plt.xlabel('number of OpenMP threads')
    plt.grid(which='both', linestyle='dashed')
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.4))
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig('timing_speedup_plot.' + fmt, bbox_inches='tight')
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
