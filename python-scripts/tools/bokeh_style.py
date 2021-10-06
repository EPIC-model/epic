from bokeh.palettes import *

bokeh_style = {
    "font.size": "25pt",
    "font.style": "normal",
    "font.font": "helvetica",
    "formatter.use_scientific": False,
    "formatter.power_limit_low": 0,
    "formatter.power_limit_high": 0,
    "formatter.precision": "auto",
}


def _reverse_palette(palette):
    return tuple(reversed(list(palette)))


Cividis256_r = _reverse_palette(Cividis256)
Inferno256_r = _reverse_palette(Inferno256)
Magma256_r = _reverse_palette(Magma256)
Plasma256_r = _reverse_palette(Plasma256)
Viridis256_r = _reverse_palette(Viridis256)

bokeh_palettes = {
    "cividis": Cividis256,
    "cividis_r": Cividis256_r,
    "inferno": Inferno256,
    "inferno_r": Inferno256_r,
    "magma": Magma256,
    "magma_r": Magma256_r,
    "plasma": Plasma256,
    "plasma_r": Plasma256_r,
    "viridis": Viridis256,
    "viridis_r": Viridis256_r,
}

# matplotlib colormaps
# 5 October 2021
# https://stackoverflow.com/questions/35315259/using-colormap-with-bokeh-scatter
# https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib
# https://matplotlib.org/stable/api/colors_api.html
# https://stackoverflow.com/questions/17873384/how-to-deep-copy-a-list
from matplotlib import cm, colors
import copy
import numpy as np

c = np.linspace(0, 1, 256)
for key in [
    "Greys",
    "Purples",
    "Blues",
    "Greens",
    "Oranges",
    "Reds",
    "YlOrBr",
    "YlOrRd",
    "OrRd",
    "PuRd",
    "RdPu",
    "BuPu",
    "GnBu",
    "PuBu",
    "YlGnBu",
    "PuBuGn",
    "BuGn",
    "YlGn",
]:
    hex_list = [None] * 256
    cmap = cm.get_cmap(key)
    for i in range(256):
        hex_list[i] = colors.to_hex(cmap(c[i]))
    bokeh_palettes["mpl_" + key] = copy.deepcopy(hex_list)

for key in [
    "PiYG",
    "PRGn",
    "BrBG",
    "PuOr",
    "RdGy",
    "RdBu",
    "RdYlBu",
    "RdYlGn",
    "Spectral",
    "coolwarm",
    "bwr",
    "seismic",
]:
    hex_list = [None] * 256
    cmap = cm.get_cmap(key)
    for i in range(256):
        hex_list[i] = colors.to_hex(cmap(c[i]))
    bokeh_palettes["mpl_" + key] = copy.deepcopy(hex_list)


# add own colormaps
from tools.colors import epic_cmaps

for key in epic_cmaps.keys():
    hex_list = [None] * 256
    cmap = epic_cmaps[key]
    for i in range(256):
        hex_list[i] = colors.to_hex(cmap(c[i]))
    bokeh_palettes[key] = copy.deepcopy(hex_list)


try:
    import colorcet as cc

    # 3 October 2021
    # https://github.com/holoviz/colorcet
    # https://colorcet.holoviz.org/getting_started/index.html
    bokeh_palettes["fire"] = cc.fire
    bokeh_palettes["bgy"] = cc.bgy
    bokeh_palettes["bgyw"] = cc.bgyw
    bokeh_palettes["bmy"] = cc.bmy
    bokeh_palettes["gray"] = cc.gray
    bokeh_palettes["kbc"] = cc.kbc
    bokeh_palettes["coolwarm"] = cc.coolwarm
    bokeh_palettes["blues"] = cc.blues
    bokeh_palettes["kb"] = cc.kb
    bokeh_palettes["kg"] = cc.kg
    bokeh_palettes["kr"] = cc.kr
    bokeh_palettes["rainbow"] = cc.rainbow
    bokeh_palettes["isolum"] = cc.isolum
    bokeh_palettes["colorwheel"] = cc.colorwheel

except:
    pass
