import matplotlib.pyplot as plt
import colorcet as cc
import matplotlib as mpl
import numpy as np
import os
import matplotlib.colors as mpl_colors
from tools.units import *

units['time'] = None

# 28 July 2022
# https://stackoverflow.com/questions/42086276/get-default-line-colour-cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

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
                cmap=cc.cm['rainbow4'], colorbar=True, cmap_norm=None, **kwargs):

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

    vmax = kwargs.pop('vmax', None)
    vmin = kwargs.pop('vmin', None)
    
    if cmap_norm == 'centered':
        norm = mpl_colors.CenteredNorm(vcenter=0.0)
    elif cmap_norm == 'symlog':
        norm = mpl_colors.SymLogNorm(linthresh=1, base=10)
    elif cmap_norm == 'log':
        norm = mpl_colors.LogNorm(vmin=vmin, vmax=vmax)
        vmin=None
        vmax=None
    else:
        norm = None

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
