import matplotlib as mpl

mpl.rcParams['figure.figsize'] = 9, 6
mpl.rcParams['figure.dpi'] = 200

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 20

legend_dict = {
    'bbox_to_anchor': (0.5, 1.3),
    'ncol':           4,
    'loc':            'upper center'
}

bokeh_style = {
    'font.size'                     : '25pt',
    'font.style'                    : 'normal',
    'font.font'                     : 'helvetica',
    'formatter.use_scientific'      : False,
    'formatter.power_limit_low'     : 0,
    'formatter.power_limit_high'    : 0,
    'formatter.precision'           : 'auto'
}
