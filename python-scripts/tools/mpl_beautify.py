import matplotlib as mpl
from tools.mpl_style import *
from tools.units import *


def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.3f"):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

    label = r"t = " + fmt % (time)
    if units['time'] is not None:
        label = r"t = \SI{" + fmt % (time) + r"}{" + units["time"] + r"}"

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox
    )


def add_number_of_parcels(plt, num, xy=(0.01, 1.05)):
    bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    plt.annotate(
        # 4 Feb 2022
        # https://tex.stackexchange.com/questions/375888/to-make-the-numero-sign-would-n-textsuperscript-underline-scriptsize-o-be-a/613930#613930
        r'N\textsuperscript{\underline{\scriptsize o}} parcels = %7.0f' % (num),
        xy=xy,
        xycoords="axes fraction",
        bbox=bbox,
    )


def add_box(plt, label, value, unit="", xy=(0.01, 1.05), fmt="%1.3f"):
    bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)


    text = label + " = " + fmt % (value) + unit
    if not unit == "":
        text = label + " = \SI{" + fmt % (value) + r"}{" + unit + r"}"

    plt.annotate(
        text, xy=xy, xycoords="axes fraction", bbox=bbox
    )


# 25 June 2021
# https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_and_donut_labels.html#sphx-glr-gallery-pie-and-polar-charts-pie-and-donut-labels-py
def get_autopct(pct, allvals):
    import numpy as np

    absolute = pct / 100.0 * np.sum(allvals)
    return "{:.2f}$\%$\n({:.1f} s)".format(pct, absolute)
