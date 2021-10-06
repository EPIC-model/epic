from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
import numpy as np

# 6 October 2021
# https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
# https://stackoverflow.com/questions/11647261/create-a-colormap-with-white-centered-around-zero
# https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.TwoSlopeNorm.html#matplotlib.colors.TwoSlopeNorm
epic_orbl = LinearSegmentedColormap.from_list(name='EpicOrangeBlue',
                                              # RGBA (red, green, blue, alpha)
                                              colors =[(1, 0.549019608, 0, 1), # darkorange
                                                       (1, 1., 1, 1),
                                                       (0.325490196, 0.2, 0.929411765, 1)], # royalblue
                                              N=255
                                              )


tol_precip_colors = [
    "#90C987",
    "#4EB256",
    "#7BAFDE",
    "#6195CF",
    "#F7CB45",
    "#EE8026",
    "#DC050C",
    "#A5170E",
    "#72190E",
    "#882E72",
    "#000000"
]

precip_colormap = ListedColormap(tol_precip_colors)

# own colormaps
epic_cmaps = {
    'epic_orbl' : epic_orbl,
    'epic_precip' : precip_colormap
}
