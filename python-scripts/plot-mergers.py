#!/usr/bin/env python
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import argparse, os
from tools.plot_style import *
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = 6, 6
mpl.rcParams['font.size'] = 16

def get_step_string(step):
    return 'step#' + str(step).zfill(10)

def get_merger_string(num):
    return 'merger#' + str(num).zfill(10)

def get_angles(a2, B11, B12, B22):
    evec = np.array([a2 - B22, B12])

    for i in range(evec.shape[1]):
        if abs(evec[0, i]) + abs(evec[1, i]) == 0.0:
            if B11[i] > B22[i]:
                evec[0, i] = evec[0, i] + np.finfo(np.float64).eps
            else:
                evec[1, i] = evec[1, i] + np.finfo(np.float64).eps

    evec = evec / np.linalg.norm(evec, 2)
    return np.arctan2(evec[1, :], evec[0, :])


def get_ellipses(pos, V, B11, B12):
    B22 = ((V / np.pi) ** 2 + B12 ** 2) / B11
    a2 = 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)

    evec = np.array([a2 - B22, B12])

    angle = get_angles(a2, B11, B12, B22)

    b2 = (V / np.pi) ** 2 / a2
    return [Ellipse(xy=pos[:, i],
                    width=2 * np.sqrt(a2[i]),
                    height=2 * np.sqrt(b2[i]),
                    angle=np.rad2deg(angle[i]))
            for i in range(len(V))]


def get_angle(a2, B11, B12, B22):
    evec = np.array([a2 - B22, B12])

    if abs(evec[0]) + abs(evec[1]) == 0.0:
        if B11[i] > B22[i]:
            evec[0] = evec[0] + np.finfo(np.float64).eps
        else:
            evec[1] = evec[1] + np.finfo(np.float64).eps

    evec = evec / np.linalg.norm(evec, 2)
    return np.arctan2(evec[1], evec[0])

def get_ellipse(pos, V, B11, B12):
    B22 = ((V / np.pi) ** 2 + B12 ** 2) / B11
    a2 = 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)

    evec = np.array([a2 - B22, B12])

    angle = get_angle(a2, B11, B12, B22)

    b2 = (V / np.pi) ** 2 / a2
    return [Ellipse(xy=pos,
                    width=2 * np.sqrt(a2),
                    height=2 * np.sqrt(b2),
                    angle=np.rad2deg(angle))]

def get_position(lower, k, h):
    return lower + k * h

def show_frames(h5merger, h5mergee, h5neighbours):

    def draw(ax, lam=None):
        global step
        #xlim = ax.get_xlim()
        #ylim = ax.get_ylim()
        ax.clear()

        s = get_step_string(step)

        # mergers
        Bm = np.array(h5merger[s]['B'])[:, nmerge]
        Vm = np.array(h5merger[s]['volume'])[nmerge]
        posm = np.array(h5merger[s]['position'])[:, nmerge]
        ell_merger = get_ellipse(posm, Vm, Bm[0], Bm[1])

        center = ell_merger[0].get_center()

        # figure out cell index (ix, iz)
        ix, iz = np.floor((center - origin) / dx)

        ix = int(ix)
        iz = int(iz)


        ixlo = ix - 1
        ixhi = ix + 2

        if ixlo < 0:
            ixlo = ix

        if ixhi > ncells[0] - 1:
            ixhi = ix

        izlo = iz - 1
        izhi = iz + 2

        if izlo < 0:
            izlo = 0

        if izhi > ncells[1]:
            izhi = ncells[1]

        xlo = get_position(origin[0], ixlo, dx[0])
        xhi = get_position(origin[0], ixhi, dx[0])
        xlim = [xlo, xhi]

        zlo = get_position(origin[1], izlo, dx[1])
        zhi = get_position(origin[1], izhi, dx[1])
        ylim = [zlo, zhi]

        ax.axhline(get_position(origin[1], iz, dx[1]),
                   color='black', linestyle='dashed', linewidth=0.75)
        ax.axhline(get_position(origin[1], iz + 1, dx[1]),
                   color='black', linestyle='dashed', linewidth=0.75)

        ax.axvline(get_position(origin[0], ix, dx[0]),
                   color='black', linestyle='dashed', linewidth=0.75)
        ax.axvline(get_position(origin[0], ix + 1, dx[0]),
                   color='black', linestyle='dashed', linewidth=0.75)

        for j, e in enumerate(ell_merger):
            ax.add_artist(e)
            e.set_alpha(0.75)
            e.set_facecolor('blue')

        # mergees
        sn = get_merger_string(nmerge)
        Bm = np.array(h5mergee[s][sn]['B'])
        Vm = np.array(h5mergee[s][sn]['volume'])
        posm = np.array(h5mergee[s][sn]['position'])

        # periodic shift
        for i in range(len(Vm)):
            if abs(posm[0, i] - center[0]) > 0.5 * extent[0]:
                if (posm[0, i] > center[0]):
                    posm[0, i] -= extent[0]
                else:
                    posm[0, i] += extent[0]

        ell_mergee = None
        ell_mergee = get_ellipses(posm, Vm, Bm[0, :], Bm[1, :])

        for j, e in enumerate(ell_mergee):
            ax.add_artist(e)
            e.set_alpha(0.75)
            e.set_facecolor('red')

        # 7 August 2021
        # https://stackoverflow.com/questions/32630818/label-plotted-ellipses
        ax.legend([ell_merger[0], ell_mergee[0]], [r'mergers', r'mergees'],
                loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.2))

        # neighbours
        Bm = np.array(h5neighbours[s][sn]['B'])
        Vm = np.array(h5neighbours[s][sn]['volume'])
        posm = np.array(h5neighbours[s][sn]['position'])

        # periodic shift
        for i in range(len(Vm)):
            if abs(posm[0, i] - center[0]) > 0.5 * extent[0]:
                if (posm[0, i] > center[0]):
                    posm[0, i] -= extent[0]
                else:
                    posm[0, i] += extent[0]

        ell_neighbours = get_ellipses(posm, Vm, Bm[0, :], Bm[1, :])

        for j, e in enumerate(ell_neighbours):
            ax.add_artist(e)
            e.set_alpha(0.5)
            e.set_facecolor('none')
            e.set_edgecolor('gray')
            e.set_linewidth(0.5)


        #
        # all other mergers and mergees
        #
        for nn in range(nmergers):
            if nn == nmerge:
                continue

            sn = get_merger_string(nn)
            Bm = np.array(h5mergee[s][sn]['B'])
            Vm = np.array(h5mergee[s][sn]['volume'])
            posm = np.array(h5mergee[s][sn]['position'])

            # periodic shift
            for i in range(len(Vm)):
                if abs(posm[0, i] - center[0]) > 0.5 * extent[0]:
                    if (posm[0, i] > center[0]):
                        posm[0, i] -= extent[0]
                    else:
                        posm[0, i] += extent[0]

            ell_mergee = None
            ell_mergee = get_ellipses(posm, Vm, Bm[0, :], Bm[1, :])

            for j, e in enumerate(ell_mergee):
                ax.add_artist(e)
                e.set_alpha(0.5)
                e.set_facecolor('gray')

        ax.set_aspect('equal', 'box')

        #bbox = dict(boxstyle="round", facecolor='wheat', alpha=0.5)
        #ax.annotate('step %3i'%(step) + ' merge %3i'%(nmerge+1) + ' out of %3i'%nmergers,
                    #xy=(0.0, 1.05),
                    #xycoords='axes fraction',
                    #bbox=bbox)
        print ("Step", step, "merge", nmerge + 1, "out of", nmergers)

        if not lam is None:
            ax.annotate(r'$\lambda =$ %1.3f'%(lam),
                    xy=(0.8, 1.05),
                    xycoords='axes fraction',
                    bbox=bbox)

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        #ax.set_xlim([origin[0], origin[0] + extent[0]])
        #ax.set_ylim([origin[1], origin[1] + extent[1]])
        ax.axvline(origin[0], linestyle='dashed', color='black')
        ax.axvline(origin[0] + extent[0], linestyle='dashed', color='black')
        ax.axhline(origin[1], linestyle='dashed', color='black')
        ax.axhline(origin[1] + extent[1], linestyle='dashed', color='black')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.tight_layout()

    def draw_frame(event, lam=None, save=False):
        ax = event.canvas.figure.get_axes()[0]
        draw(ax, lam)
        #if save:
            #fname = '/home/matthias/Documents/projects/papers/epic2d/data/Straka/lambda_mergee/'
            #fname = fname + 'step_' + str(step) + '_merge_' + str(nmerge) + '.png'
            #fig.canvas.print_png(fname)
            #exit()
        #else:
        fig.canvas.draw()
        return

    def on_press(event):
        global keys, step, nmerge, nmergers
        if event.key not in '':
            keys.add(event.key)
        if len(keys) == 2 and '+' in keys:
            nmerge = 0
            step = step + 1
            try:
                s = get_step_string(step)
                nmergers = h5mergee[s].attrs['nmergers'][0]
                Vm = np.array(h5merger[s]['volume'])
                if len(Vm) != nmergers:
                    print("Number of mergers", nmergers, "!=", len(Vm))
                    exit()
            except:
                print("Beyond maximum step number. Start from the beginning!")
                step = 0

            draw_frame(event)
        if len(keys) == 1 and 'c' in keys:
            found = False
            nmerge = nmerge + 1
            while found == False:
                s = get_step_string(step)
                nmergers = h5mergee[s].attrs['nmergers'][0]
                for nnn in range(nmerge, nmergers):
                    print("Merge", nnn)
                    sn = get_merger_string(nnn)
                    Vm = np.array(h5mergee[s][sn]['volume'])
                    if Vm.shape[0] > 3:
                        print(nnn, Vm.shape[0])
                        found = True
                        nmerge = nnn
                        break
                    if Vm.shape[0] == '':
                        exit()
                if found == False:
                    step = step + 1
                    nmerge = 0
                    print("Step", step)

            draw_frame(event)

        #if len(keys) == 1 and 'l' in keys:
            #found = False
            #nmerge = nmerge + 1
            #while found == False:
                #s = get_step_string(step)
                #nmergers = h5mergee[s].attrs['nmergers'][0]
                #for nnn in range(nmerge, nmergers):
                    #print("Merge", nnn)
                    #Bm = np.array(h5merger[s]['B'])[:, nnn]
                    #Vm = np.array(h5merger[s]['volume'])[nnn]

                    #B11 = Bm[0]
                    #B12 = Bm[1]
                    #B22 = ((Vm / np.pi) ** 2 + B12 ** 2) / B11
                    #a2 = 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)
                    #ab = Vm / np.pi
                    #lam = a2 / ab

                    #if lam > 4:
                        #print(nnn, lam)
                        #found = True
                        #nmerge = nnn
                        #break
                #if found == False:
                    #step = step + 1
                    #nmerge = 0
                    #print("Step", step)
            #print("Found", step, nmerge, lam)
            #draw_frame(event, lam)

        #if len(keys) == 1 and 'm' in keys:
            #found = False
            #nmerge = nmerge + 1
            #while found == False:
                #s = get_step_string(step)
                #nmergers = h5mergee[s].attrs['nmergers'][0]
                #for nnn in range(nmerge, nmergers):
                    #print("Merge", nnn)
                    #sn = get_merger_string(nnn)
                    #Bm = np.array(h5mergee[s][sn]['B'])
                    #Vm = np.array(h5mergee[s][sn]['volume'])
                    #lam = 0
                    #for ii in range(Vm.shape[0]):
                        #B11 = Bm[0, ii]
                        #B12 = Bm[1, ii]
                        #B22 = ((Vm[ii] / np.pi) ** 2 + B12 ** 2) / B11
                        #a2 = 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)
                        #ab = Vm[ii] / np.pi
                        #lam = max(lam, a2 / ab)

                    #if lam > 4:
                        #print(nnn, lam)
                        #found = True
                        #nmerge = nnn
                        #break
                #if found == False:
                    #step = step + 1
                    #nmerge = 0
                    #print("Step", step)
            #print("Found", step, nmerge, lam)
            #draw_frame(event, lam, save=True)

        if len(keys) == 1 and 'n' in keys:
            nmerge = nmerge + 1
            if nmerge == nmergers:
                step = step + 1
                nmerge = 0
                s = get_step_string(step)
                nmergers = h5mergee[s].attrs['nmergers'][0]
                Vm = np.array(h5merger[s]['volume'])
                if len(Vm) != nmergers:
                    print("Number of mergers", nmergers, "!=", len(Vm))
                    exit()

            draw_frame(event)
        if len(keys) == 1 and 'p' in keys:
            nmerge = nmerge - 1

            if nmerge < 0:
                nmerge = 0
                step = step - 1
                if step < 0:
                    step = 0
                    print("Negative step number. Go to step 0.")
                s = get_step_string(step)
                nmergers = h5mergee[s].attrs['nmergers'][0]
                nmerge = nmergers - 1
                Vm = np.array(h5merger[s]['volume'])
                if len(Vm) != nmergers:
                    print("Number of mergers", nmergers, "!=", len(Vm))
                    exit()
            draw_frame(event)
        if len(keys) == 1 and '-' in keys:
            step = step - 1
            if step == -1:
                print("Negative step number. Go to step 0.")
                step = 0
            s = get_step_string(step)
            nmergers = h5mergee[s].attrs['nmergers'][0]
            nmerge = nmergers - 1
            Vm = np.array(h5merger[s]['volume'])
            if len(Vm) != nmergers:
                print("Number of mergers", nmergers, "!=", len(Vm))
                exit()

            draw_frame(event)
        if 'r' in keys:
            keys.clear()
            keys.add('r')
        return

    def on_release(event):
        global keys
        keys.remove(event.key)
        return

    global origin, extent, step, keys, nmergers, nmerge, dx
    step = 0
    nmerge = 0

    s = get_step_string(step)
    nmergers = h5mergee[s].attrs['nmergers'][0]
    origin = np.array(h5merger['box'].attrs['origin'])
    extent = np.array(h5merger['box'].attrs['extent'])
    ncells = np.array(h5merger['box'].attrs['ncells'])

    dx = extent / ncells

    keys=set()

    fig = plt.figure(num=1, figsize=(6, 6), dpi=200)
    ax = plt.gca()

    draw(ax)

    cid=fig.canvas.mpl_connect('key_press_event',on_press)
    cid2=fig.canvas.mpl_connect('key_release_event',on_release)

    plt.show()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Plot individual mergers and their mergees. " + \
            "EPIC must be compiled with '--enable-merger-dump'.")

    parser.add_argument("--mergers",
                          type=str,
                          required=True,
                          help="*_merger.hdf5 output file of EPIC")

    parser.add_argument("--mergees",
                        type=str,
                        required=True,
                        help="*_mergee.hdf5 output file of EPIC")

    parser.add_argument("--neighbours",
                        type=str,
                        required=True,
                        help="*_neighbours.hdf5 output file of EPIC")

    args = parser.parse_args()

    if not os.path.exists(args.mergers):
        raise IOError("File '" + args.mergers + "' does not exist.")

    if not os.path.exists(args.mergees):
        raise IOError("File '" + args.mergees + "' does not exist.")

    if not os.path.exists(args.neighbours):
        raise IOError("File '" + args.neighbours + "' does not exist.")

    print("--------------------------------------------------------")
    print("Navigation:")
    print("\t +     move to next time step")
    print("\t -     move to previous time step")
    print("\t n     move to next merge within time step")
    print("\t p     move to prevous merge within time step")
    print("--------------------------------------------------------")

    fmerger = args.mergers
    fmergee = args.mergees
    fneighbours = args.neighbours

    h5merger = h5py.File(fmerger, 'r')
    h5mergee = h5py.File(fmergee, 'r')
    h5neighbours = h5py.File(fneighbours, 'r')

    step = 0

    show_frames(h5merger, h5mergee, h5neighbours)

    h5merger.close()
    h5mergee.close()
    h5neighbours.close()
