from bokeh.palettes import *

bokeh_style = {
    'font.size'                     : '25pt',
    'font.style'                    : 'normal',
    'font.font'                     : 'helvetica',
    'formatter.use_scientific'      : False,
    'formatter.power_limit_low'     : 0,
    'formatter.power_limit_high'    : 0,
    'formatter.precision'           : 'auto'
}


def _reverse_palette(palette):
    return tuple(reversed(list(palette)))


Cividis256_r = _reverse_palette(Cividis256)
Inferno256_r = _reverse_palette(Inferno256)
Magma256_r = _reverse_palette(Magma256)
Plasma256_r  = _reverse_palette(Plasma256)
Viridis256_r = _reverse_palette(Viridis256)

bokeh_palettes = {
    'cividis'           : Cividis256,
    'cividis_r'         : Cividis256_r,
    'inferno'           : Inferno256,
    'inferno_r'         : Inferno256_r,
    'magma'             : Magma256,
    'magma_r'           : Magma256_r,
    'plasma'            : Plasma256,
    'plasma_r'          : Plasma256_r,
    'viridis'           : Viridis256,
    'viridis_r'         : Viridis256_r
}


try:
    import colorcet as cc

    # 3 October 2021
    # https://github.com/holoviz/colorcet
    # https://colorcet.holoviz.org/getting_started/index.html
    bokeh_palettes['fire'] = cc.fire
    bokeh_palettes['bgy'] = cc.bgy
    bokeh_palettes['bgyw'] = cc.bgyw
    bokeh_palettes['bmy'] = cc.bmy
    bokeh_palettes['gray'] = cc.gray
    bokeh_palettes['kbc'] = cc.kbc
    bokeh_palettes['coolwarm'] = cc.coolwarm
    bokeh_palettes['blues'] = cc.blues
    bokeh_palettes['kb'] = cc.kb
    bokeh_palettes['kg'] = cc.kg
    bokeh_palettes['kr'] = cc.kr
    bokeh_palettes['rainbow'] = cc.rainbow
    bokeh_palettes['isolum'] = cc.isolum
    bokeh_palettes['colorwheel'] = cc.colorwheel

except:
    pass
