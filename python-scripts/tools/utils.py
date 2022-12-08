import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
import numpy as np
import os
import matplotlib.colors as mpl_colors

mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 11,
    "text.usetex": True,
    'legend.framealpha': 1.0,
    'lines.linewidth': 0.75,
    'grid.color':     'b0b0b0',
    'grid.linestyle': 'dotted',
    'grid.linewidth': 0.25,
    'grid.alpha':     1.0,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{bm}"
#        r"\usepackage{siunitx}",
        ])
})

# 28 July 2022
# https://stackoverflow.com/questions/42086276/get-default-line-colour-cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.2f", **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    label = r"$t = " + fmt % (time) + "$"
    plt.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox, **kwargs)

def add_annotation(ax, label, xy, **kwargs):
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
    ax.annotate(label, xy=xy, xycoords="axes fraction", bbox=bbox, **kwargs)

def remove_xticks(ax):
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

def remove_yticks(ax):
    ax.tick_params(axis='y', which='both', right=False, left=False)

# assumes data in ordering (x, y, z)
def get_max_magnitude(xcomp, ycomp, zcomp, plane, loc):
    if plane == 'xy':
        magn = xcomp[:, :, loc] ** 2 + ycomp[:, :, loc] ** 2 + zcomp[:, :, loc] ** 2
    elif plane == 'xz':
        magn = xcomp[:, loc, :] ** 2 + ycomp[:, loc, :] ** 2 + zcomp[:, loc, :] ** 2
    elif plane == 'yz':
        magn = xcomp[loc, :, :] ** 2 + ycomp[loc, :, :] ** 2 + zcomp[loc, :, :] ** 2
    return np.sqrt(magn.max())

def get_plane(plane, loc, fdata):
    if plane == 'xy':
        pl = fdata[:, :, loc]
    elif plane == 'xz':
        pl = fdata[:, loc, :]
    elif plane == 'yz':
        pl = fdata[loc, :, :]
    return pl

# assumes fdata ordering (x, y, z)
def make_imshow(ax, plane, loc, fdata, ncr,
                cmap='rainbow4', colorbar=True, cmap_norm=None, **kwargs):

    if ncr == None:
        origin = kwargs.get('origin')
        extent = kwargs.get('extent')
    else:
        origin = ncr.get_box_origin()
        extent = ncr.get_box_extent()

    # plane = 'xy':
    imin = origin[0]
    imax = origin[0] + extent[0]
    jmin = origin[1]
    jmax = origin[1] + extent[1]

    if plane == 'xy':
        pl = fdata[:, :, loc]
        xlab = r'$x$'
        ylab = r'$y$'
    elif plane == 'xz':
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = fdata[:, loc, :]
        xlab = r'$x$'
        ylab = r'$z$'
    elif plane == 'yz':
        imin = jmin
        imax = jmax
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = fdata[loc, :, :]
        xlab = r'$y$'
        ylab = r'$z$'

    if cmap_norm == 'centered':
        norm = mpl_colors.CenteredNorm(vcenter=0.0)
    elif cmap_norm == 'symlog':
        norm = mpl_colors.SymLogNorm(linthresh=1, base=10)
    elif cmap_norm == 'log':
        norm = mpl_colors.LogNorm()
    else:
        norm = None

        
    vmax = kwargs.pop('vmax', None)
    vmin = kwargs.pop('vmin', None)

    cmap = cc.cm[cmap]
    color_under = kwargs.pop('color_under', 'darkblue')
    color_over = kwargs.pop('color_over', 'orangered')

    cbar_ext = kwargs.pop('cbar_ext', False)

    if cbar_ext:
        cmap.set_under(color=color_under)
        cmap.set_over(color=color_over)

    im = ax.imshow(X=pl.transpose(),
                   vmin=vmin,
                   vmax=vmax,
                   cmap=cmap,
                   norm=norm,
                   interpolation='bilinear',
                   origin='lower',
                   extent=[imin, imax, jmin, jmax])

    xticks = kwargs.pop('xticks', None)
    xticklab = kwargs.pop('xticklab', None)

    if not xticks is None and not xticklab is None:
        ax.set_xticks(xticks, xticklab)

    yticks = kwargs.pop('yticks', None)
    yticklab = kwargs.pop('yticklab', None)
    
    if not yticks is None and not yticklab is None:
        ax.set_yticks(yticks, yticklab)
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    cbar = None
    if colorbar:
        if cbar_ext:
            cbar = ax.cax.colorbar(im, extend='both', extendrect=False)
        else:
            cbar = ax.cax.colorbar(im)

        clabel = kwargs.pop('clabel', None)
        if not clabel is None:
            cbar.set_label(clabel)
            
        if not cmap_norm == 'symlog' and not cmap_norm == 'log':
            cbar.formatter.set_powerlimits((0, 0))
            # 18 July 2022
            # https://stackoverflow.com/questions/34039396/matplotlib-colorbar-scientific-notation-offset
            cbar.ax.yaxis.set_offset_position('left')
    return im, cbar

def make_mean_profiles(ax, ncr, step, fields, labels, normalise=False, **kwargs):
    z = ncr.get_all('z')
    n = len(z)
    zticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    zticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

    xticks = kwargs.pop('xticks', [-1, 0, 1])
    xlim = kwargs.pop('xlim', [-1.1, 1.1])

    for i, field in enumerate(fields):
        bar = np.zeros(n)
        data = ncr.get_dataset(step=step, name=field, copy_periodic=False)
        bar = data.mean(axis=(0, 1))

        if normalise:
            bar = bar / bar.max()

        ax.plot(bar, z,
                color=colors[i],
                label=labels[i])
        ax.grid(zorder=-1)
        #ax.set_aspect(1)
        ax.set_yticks(ticks=zticks, labels=zticklab)
        #ax.set_xticks(ticks=xticks)
        #ax.set_xlim(xlim)
    return ax

def make_rms_profiles(ax, ncr, step, fields, labels):
    z = ncr.get_all('z')
    n = len(z)
    zticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    zticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

    for i, field in enumerate(fields):
        rms = np.zeros(n)
        data = ncr.get_dataset(step=step, name=field, copy_periodic=False)
        rms = np.sqrt((data ** 2).mean(axis=(0, 1)))

        ax.plot(rms, z,
                color=colors[i],
                label=labels[i])
        ax.grid(zorder=-1)
        ax.set_yticks(ticks=zticks, labels=zticklab)
    return ax


def save_figure(plt, figpath, fignum=1, overwrite=False):
    figname = 'fig' + str(fignum) + '.pdf'
    fname = os.path.join(figpath, figname)

    if os.path.exists(fname) and not overwrite:
        print("Figure '" + fname + "' already exists.")
        plt.close()
        exit()

    print("Save figure as:", fname)
    plt.savefig(fname=fname, format='pdf', bbox_inches='tight')
    plt.close()

# 7 July 2022
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(t, tref):
    idx = np.argmin(np.abs(tref - t))
    return idx
