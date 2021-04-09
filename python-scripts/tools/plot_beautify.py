import matplotlib as mpl

mpl.rcParams['text.usetex'] = True

def add_timestamp(plt, time):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor='wheat', alpha=0.5)
    plt.annotate('t = %.3f'%(time),
                 xy=(0.78, 1.05),
                 xycoords='axes fraction',
                 bbox=bbox)

def add_number_of_parcels(plt, num):
    bbox = dict(boxstyle="round", facecolor='wheat', alpha=0.5)
    plt.annotate('no. ellipses = %7.0f'%(num),
                 xy=(0.01, 1.05),
                 xycoords='axes fraction',
                 bbox=bbox)

def add_box(plt, label, value, unit='', xy=(0.01, 1.05), fmt='%1.3f'):
    bbox = dict(boxstyle="round", facecolor='wheat', alpha=0.5)
    plt.annotate(label + ' = ' + fmt%(value) + unit,
                 xy=xy,
                 xycoords='axes fraction',
                 bbox=bbox)
