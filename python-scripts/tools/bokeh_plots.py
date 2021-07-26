from bokeh.io import export_png, export_svg
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, ColorBar
from bokeh.palettes import Viridis256
from bokeh.transform import linear_cmap
from tools.h5_reader import H5Reader
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

    font_size = '15pt'

    if display == 'full HD':
        pw = 1920
        ph = 1080
    elif display == 'UHD':
        font_size = '25pt'
        pw = 3840
        ph = 2160
    elif display == '4K':
        font_size = '32pt'
        pw = 4096
        ph = 2160
    else:
        pw = np.nanmin([1920,int(1080*(right-left)/(top-bottom))])
        ph = np.nanmin([1080,int(1920*(top-bottom)/(right-left))])

    nparcels= h5reader.get_num_parcels(step)
    ttime = h5reader.get_step_attribute(step=step, name='t')

    graph = figure(output_backend="webgl",
                   plot_width=pw,
                   plot_height = ph,
                   aspect_ratio = (right-left)/(top-bottom),
                   x_range = (left, right),
                   y_range = (bottom, top),
                   x_axis_label='x',
                   y_axis_label='z',
                   title = coloring + \
                       '                              time = %15.3f'%ttime + \
                       '                              #parcels = %10d'%nparcels)

    # 20 July 2021
    # https://stackoverflow.com/questions/32158939/python-bokeh-remove-toolbar-from-chart
    graph.toolbar.logo = None
    graph.toolbar_location = None

    # 20 July 2021
    # https://stackoverflow.com/questions/47220491/how-do-you-change-ticks-label-sizes-using-pythons-bokeh
    graph.xaxis.axis_label_text_font_size = font_size
    graph.xaxis.major_label_text_font_size = font_size
    graph.yaxis.axis_label_text_font_size = font_size
    graph.yaxis.major_label_text_font_size = font_size

    graph.title.text_font_size = font_size

    x, y, width, height, angle = h5reader.get_ellipses_for_bokeh(step)

    if coloring == 'aspect-ratio':
        data = h5reader.get_aspect_ratio(step=step)
    else:
        data = h5reader.get_dataset(step=step, name=coloring)


    source = ColumnDataSource(dict(x=x,y=y, width=width, height=height,
                                   angle = angle, fill_color=data))

    Viridis256_r = tuple(reversed(list(Viridis256)))
    mapper = linear_cmap(field_name='fill_color', palette=Viridis256_r,
                         low=vmin, high=vmax)

    graph.ellipse(x='x', y='y', width='width', height='height',angle='angle',
    color = mapper,fill_alpha=0.75,line_color=None,source=source)
    color_bar = ColorBar(color_mapper=mapper['transform'], label_standoff=12,
                         title_text_font_size=font_size,
                         major_label_text_font_size=font_size)
    graph.add_layout(color_bar, 'right')

    return graph


def bokeh_plot_parcels(fname, step, shw=False, fmt='png',
                       coloring='aspect-ratio', display='full HD'):

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
    else:
        vmin, vmax = h5reader.get_dataset_min_max(coloring)

    graph = _bokeh_plot_parcels(h5reader, step, coloring, vmin, vmax, display)

    if shw:
        show(graph)
    elif fmt == 'png':
        export_png(graph, filename = 'parcels_'  + coloring + '_step_' + \
            str(step).zfill(len(str(nsteps))) + '.png')
    elif fmt == 'svg':
        export_svg(graph, filename = 'parcels_'  + coloring + '_step_' + \
            str(step).zfill(len(str(nsteps))) + '.svg')
    else:
        raise IOError("Bokeh plot does not support '" + fmt + "' format.")

    h5reader.close()
