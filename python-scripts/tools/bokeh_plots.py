from bokeh.io import export_png, export_svg
import bokeh.plotting as bpl
from bokeh.models import ColumnDataSource, ColorBar, FixedTicker

from bokeh.palettes import Viridis256
from bokeh.transform import linear_cmap
from tools.h5_reader import H5Reader
from tools.plot_style import *
import numpy as np

def _bokeh_plot_parcels(h5reader, step, coloring, vmin, vmax, display=None, **kwargs):
    nsteps = h5reader.get_num_steps()
    extent = h5reader.get_box_extent()
    origin = h5reader.get_box_origin()

    left = origin[0]
    right = origin[0] + extent[0]
    bottom = origin[1]
    top = origin[1] + extent[1]

    # instantiating the figure object
    fkwargs = {k: v for k, v in kwargs.items() if v is not None}
    left = fkwargs.get('xmin', origin[0])
    right = fkwargs.get('xmax', origin[0] + extent[0])
    bottom = fkwargs.get('ymin', origin[1])
    top = fkwargs.get('ymax', origin[1] + extent[1])

    font_size = bokeh_style['font.size']
    text_font = bokeh_style['font.font']

    if display == 'full HD':
        pw = 1920
        ph = 1080
    elif display == 'UHD':
        pw = 3840
        ph = 2160
    elif display == '4K':
        pw = 4096
        ph = 2160
    else:
        pw = np.nanmin([1920,int(1080*(right-left)/(top-bottom))])
        ph = np.nanmin([1080,int(1920*(top-bottom)/(right-left))])

    nparcels= h5reader.get_num_parcels(step)
    ttime = h5reader.get_step_attribute(step=step, name='t')

    if coloring == 'aspect-ratio':
        title = 'aspect ratio'
    elif coloring == 'vol-distr':
        title = 'volume distribution'
    else:
        title = coloring

    graph = bpl.figure(output_backend="webgl",
                       plot_width=pw,
                       plot_height = ph,
                       aspect_ratio = (right-left)/(top-bottom),
                       x_range = (left, right),
                       y_range = (bottom, top),
                       x_axis_label='x (m)',
                       y_axis_label='y (m)',
                       title = title + \
                       '\t\t\t\t time = %15.3f'%ttime + \
                       '\t\t\t\t #parcels = %10d'%nparcels)

    # 20 July 2021
    # https://stackoverflow.com/questions/32158939/python-bokeh-remove-toolbar-from-chart
    graph.toolbar.logo = None
    graph.toolbar_location = None


    # 20 July 2021
    # https://stackoverflow.com/questions/47220491/how-do-you-change-ticks-label-sizes-using-pythons-bokeh
    graph.xaxis.axis_label_text_font_size = font_size
    graph.xaxis.major_label_text_font_size = font_size
    graph.xaxis.axis_label_text_font = text_font
    graph.xaxis.major_label_text_font = text_font
    graph.yaxis.axis_label_text_font_size = font_size
    graph.yaxis.major_label_text_font_size = font_size
    graph.yaxis.axis_label_text_font = text_font
    graph.yaxis.major_label_text_font = text_font

    graph.yaxis.formatter.use_scientific = bokeh_style['formatter.use_scientific']
    graph.xaxis.formatter.use_scientific = bokeh_style['formatter.use_scientific']
    graph.yaxis.formatter.power_limit_low = bokeh_style['formatter.power_limit_low']
    graph.xaxis.formatter.power_limit_low = bokeh_style['formatter.power_limit_low']
    graph.yaxis.formatter.power_limit_high = bokeh_style['formatter.power_limit_high']
    graph.xaxis.formatter.power_limit_high = bokeh_style['formatter.power_limit_high']
    graph.yaxis.formatter.precision = bokeh_style['formatter.precision']
    graph.xaxis.formatter.precision = bokeh_style['formatter.precision']

    graph.title.text_font_size = font_size
    graph.title.text_font = text_font

    x, y, width, height, angle = h5reader.get_ellipses_for_bokeh(step)

    if coloring == 'aspect-ratio':
        data = h5reader.get_aspect_ratio(step=step)
    elif coloring == 'vol-distr':
        data = h5reader.get_dataset(step=step, name='volume')
        data[data <= vmin] = 0.0
        data[data > vmin] = 1.0
    else:
        data = h5reader.get_dataset(step=step, name=coloring)


    source = ColumnDataSource(dict(x=x,y=y, width=width, height=height,
                                   angle = angle, fill_color=data))

    mapper = None
    color_bar = None
    if coloring == 'vol-distr':
        mapper = linear_cmap(field_name='fill_color', palette=['blue', 'red'], low=0, high=1)
        # 5 August 2021
        # https://stackoverflow.com/questions/63344015/how-do-i-get-a-bokeh-colorbar-to-show-the-min-and-max-value
        ticker = FixedTicker(ticks=[0, 0.5, 1])
        color_bar = ColorBar(color_mapper=mapper['transform'], label_standoff=12,
                             ticker=ticker,
                            title_text_font_size=font_size,
                            major_label_text_font=text_font,
                            major_label_text_font_size=font_size)
        # 5 August 2021
        # https://stackoverflow.com/questions/37173230/how-do-i-use-custom-labels-for-ticks-in-bokeh
        color_bar.major_label_overrides = {0: '0', 0.5: 'Vmin', 1: 'Vmax'}
    else:
        Viridis256_r = tuple(reversed(list(Viridis256)))
        mapper = linear_cmap(field_name='fill_color', palette=Viridis256_r,
                             low=vmin, high=vmax)
        color_bar = ColorBar(color_mapper=mapper['transform'], label_standoff=12,
                            title_text_font_size=font_size,
                            major_label_text_font=text_font,
                            major_label_text_font_size=font_size)

    graph.ellipse(x='x', y='y', width='width', height='height',angle='angle',
                  color = mapper,fill_alpha=0.75,line_color=None,source=source)
    graph.add_layout(color_bar, 'right')

    return graph


def bokeh_plot_parcels(fname, step, show=False, fmt='png',
                       coloring='aspect-ratio', display='full HD', **kwargs):

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
        vmin = vcell / h5reader.get_parcel_option('min_vratio')
        vmax = vcell / h5reader.get_parcel_option('max_vratio')
    else:
        vmin, vmax = h5reader.get_dataset_min_max(coloring)

    graph = _bokeh_plot_parcels(h5reader, step, coloring, vmin, vmax, display, **kwargs)

    if show:
        bpl.show(graph)
    elif fmt == 'png':
        export_png(graph, filename = 'parcels_'  + coloring + '_step_' + \
            str(step).zfill(len(str(nsteps))) + '.png')
    elif fmt == 'svg':
        export_svg(graph, filename = 'parcels_'  + coloring + '_step_' + \
            str(step).zfill(len(str(nsteps))) + '.svg')
    else:
        raise IOError("Bokeh plot does not support '" + fmt + "' format.")

    h5reader.close()
