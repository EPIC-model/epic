from bokeh.io import export_png, export_svg
from bokeh.io.export import get_screenshot_as_png
import bokeh.plotting as bpl
from bokeh.models import ColumnDataSource, ColorBar, FixedTicker, LinearColorMapper, SingleIntervalTicker
from bokeh.transform import linear_cmap
from tools.nc_reader import nc_reader
from tools.bokeh_style import *
from tools.units import *
import numpy as np
from scipy import interpolate

def _get_bokeh_basic_graph(origin, extent, title=None, **kwargs):
    no_xaxis = kwargs.pop("no_xaxis", False)
    no_yaxis = kwargs.pop("no_yaxis", False)
    no_xlabel = kwargs.pop("no_xlabel", False)
    no_ylabel = kwargs.pop("no_ylabel", False)
    display = kwargs.pop("display", "full HD")
    plot_width = kwargs.pop('plot_width', None)
    plot_height = kwargs.pop('plot_height', None)
    yaxis_interval = kwargs.pop('yaxis_interval', None)
    xaxis_interval = kwargs.pop('xaxis_interval', None)

    # instantiating the figure object
    fkwargs = {k: v for k, v in kwargs.items() if v is not None}
    left = fkwargs.get("xmin", origin[0])
    right = fkwargs.get("xmax", origin[0] + extent[0])
    bottom = fkwargs.get("ymin", origin[1])
    top = fkwargs.get("ymax", origin[1] + extent[1])

    font_size = bokeh_style["font.size"]
    text_font = bokeh_style["font.font"]

    if plot_width is None and plot_height is None:
        if display == "full HD":
            plot_width = 1920
            plot_height = 1080
        elif display == "UHD":
            plot_width = 3840
            plot_height = 2160
        elif display == "4K":
            plot_width = 4096
            plot_height = 2160
        else:
            plot_width = np.nanmin([1920, int(1080 * (right - left) / (top - bottom))])
            plot_height = np.nanmin([1080, int(1920 * (top - bottom) / (right - left))])

    x_axis_label = get_label("x", units["position"], is_bokeh=True)
    y_axis_label = get_label("y", units["position"], is_bokeh=True)

    if no_xlabel:
        x_axis_label = " "

    if no_ylabel:
        y_axis_label = " "

    graph = bpl.figure(
        output_backend="webgl",
        plot_width=plot_width,
        plot_height=plot_height,
        aspect_ratio=(right - left) / (top - bottom),
        x_range=(left, right),
        y_range=(bottom, top),
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
        title=title,
    )

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

    if no_xaxis:
        # 27 Sept. 2021
        # https://stackoverflow.com/questions/30545372/hide-axis-in-bokeh
        # https://stackoverflow.com/questions/36112404/bokeh-gridplot-tighten-layout-and-outer-sidelabels
        # graph.xaxis.visible = False
        graph.xaxis.major_label_text_font_size = "0pt"

    if no_yaxis:
        graph.yaxis.major_label_text_font_size = "0pt"
        # graph.xaxis.visible = False

    graph.yaxis.formatter.use_scientific = bokeh_style["formatter.use_scientific"]
    graph.xaxis.formatter.use_scientific = bokeh_style["formatter.use_scientific"]
    graph.yaxis.formatter.power_limit_low = bokeh_style["formatter.power_limit_low"]
    graph.xaxis.formatter.power_limit_low = bokeh_style["formatter.power_limit_low"]
    graph.yaxis.formatter.power_limit_high = bokeh_style["formatter.power_limit_high"]
    graph.xaxis.formatter.power_limit_high = bokeh_style["formatter.power_limit_high"]
    graph.yaxis.formatter.precision = bokeh_style["formatter.precision"]
    graph.xaxis.formatter.precision = bokeh_style["formatter.precision"]

    if not yaxis_interval is None:
        graph.yaxis.ticker = SingleIntervalTicker(interval=yaxis_interval)

    if not xaxis_interval is None:
        graph.xaxis.ticker = SingleIntervalTicker(interval=xaxis_interval)

    if not title is None:
        graph.title.text_font_size = font_size
        graph.title.text_font = text_font

    return graph, (plot_width, plot_height)


def _bokeh_save(graph, fname, fmt, show, **kwargs):
    jpg_quality = kwargs.pop("jpg_quality", 90)

    if show:
        bpl.show(graph)
    elif fmt == "png":
        export_png(graph, filename=fname + ".png", timeout=120)
    elif fmt == "svg":
        export_svg(graph, filename=fname + ".svg", timeout=120)
    elif fmt == "jpg":
        # 4 October 2021
        # https://docs.bokeh.org/en/latest/docs/reference/io.html#bokeh.io.export.get_screenshot_as_png
        im = get_screenshot_as_png(graph, timeout=120)
        # 29 Sept. 2021
        # https://stackoverflow.com/questions/4353019/in-pythons-pil-how-do-i-change-the-quality-of-an-image
        # https://stackoverflow.com/questions/43258461/convert-png-to-jpeg-using-pillow
        rgb_im = im.convert("RGB")
        rgb_im.save(fname + ".jpg", quality=jpg_quality)
    else:
        raise IOError("Bokeh plot does not support '" + fmt + "' format.")


def _bokeh_plot_field(ncreader, step, field, vmin, vmax, hybrid=False, **kwargs):
    no_title = kwargs.pop("no_title", False)
    no_colorbar = kwargs.pop("no_colorbar", False)
    zoom_factor = kwargs.pop("zoom_factor", 1.0)
    norm = kwargs.pop("norm", False)
    color_bar_width = kwargs.pop('color_bar_width', 'auto')
    high_color = kwargs.pop('high_color', None)
    low_color = kwargs.pop('low_color', None)

    cmap = kwargs.get("cmap", "viridis_r")
    if not cmap in bokeh_palettes.keys():
        raise ValueError("Colormap '" + cmap + "' not available.")
    palette = bokeh_palettes[cmap]

    nsteps = ncreader.get_num_steps()
    extent = ncreader.get_box_extent()
    origin = ncreader.get_box_origin()

    ttime = ncreader.get_dataset(step=step, name="t")
    title = field + "\t\t\t\t time = %15.3f" % ttime

    if no_title:
        title = None

    graph, (pw, ph) = _get_bokeh_basic_graph(
        origin=origin, extent=extent, title=title, **kwargs
    )

    data = np.transpose(ncreader.get_dataset(step=step, name=field))

    font_size = bokeh_style["font.size"]
    text_font = bokeh_style["font.font"]



    # Shift the data to the correct position
    ny_input = np.shape(data)[0]
    nx_input = np.shape(data)[1]
    dx = extent[0] / nx_input
    dy = extent[1] / (ny_input - 1)
    data_periodic = np.zeros((ny_input, nx_input + 1))
    data_periodic[:, :nx_input] = data
    data_periodic[:, nx_input] = data[:, 0]

    # Zoom in on the data
    if hybrid:
        ncells = ncells = ncreader.get_box_ncells()
        zoom = (
            max(int(zoom_factor * pw / ncells[0]), int(zoom_factor * ph / ncells[1]))
            + 1
        )
        x_input = np.linspace(
            origin[0], origin[0] + extent[0], nx_input + 1
        )  # Include boundary
        y_input = np.linspace(origin[1], origin[1] + extent[1], ny_input)
        x_zoom = np.linspace(origin[0], origin[0] + extent[0], zoom * nx_input + 1)
        y_zoom = np.linspace(
            origin[1], origin[1] + extent[1], zoom * (ny_input - 1) + 1
        )
        interp_f = interpolate.RegularGridInterpolator(
            (x_input, y_input), data_periodic.T
        )
        # Replace the fields
        dx = x_zoom[1] - x_zoom[0]
        dy = y_zoom[1] - y_zoom[0]
        meshgrid_points = np.meshgrid(x_zoom, y_zoom)
        flat_grid = np.array([m.flatten() for m in meshgrid_points])
        out_array = interp_f(flat_grid.T)
        data_periodic = out_array.reshape(meshgrid_points[0].shape)

    ticker = None

    dmax = vmax
    dmin = vmin

    if norm:
        data_periodic[data_periodic > 0] = data_periodic[data_periodic > 0] / dmax
        data_periodic[data_periodic < 0] = data_periodic[data_periodic < 0] / (-dmin)

        vmin = -1.0
        vmax = 1.0

        ticker = FixedTicker(ticks=[-1, -0.5, 0, 0.5, 1])


    mapper = LinearColorMapper(palette=palette, low=vmin, high=vmax,
                               high_color=high_color, low_color=low_color)

    color_bar = ColorBar(
        color_mapper=mapper,
        label_standoff=12,
        title_text_font_size=font_size,
        major_label_text_font=text_font,
        major_label_text_font_size=font_size,
    )

    if not ticker is None:
        color_bar.ticker = ticker

    if norm:
        color_bar.major_label_overrides = {
            -1: str(round(dmin, 4)),
            -0.5: str(round(dmin/2, 4)),
            0: "0",
            0.5: str(round(dmax/2, 4)),
            1: str(round(dmax, 4))
        }

    graph.image(
        image=[data_periodic],
        x=origin[0] - 0.5 * dx,
        y=origin[1] - 0.5 * dy,
        dw=extent[0] + dx,
        dh=extent[1] + dy,
        color_mapper=mapper,
    )


    if not no_colorbar:
        color_bar.width = color_bar_width
        graph.add_layout(color_bar, "right")

    return graph


def _bokeh_plot_parcels(ncreader, step, coloring, vmin, vmax, **kwargs):
    no_title = kwargs.pop("no_title", False)
    no_colorbar = kwargs.pop("no_colorbar", False)
    graph = kwargs.pop("graph", None)
    fill_alpha = kwargs.pop("fill_alpha", 0.75)
    title = kwargs.pop("title", None)
    norm = kwargs.pop("norm", False)
    color_bar_width = kwargs.pop('color_bar_width', 'auto')

    cmap = kwargs.get("cmap", "viridis_r")
    if not cmap in bokeh_palettes.keys():
        raise ValueError("Colormap '" + cmap + "' not available.")
    palette = bokeh_palettes[cmap]

    nsteps = ncreader.get_num_steps()
    extent = ncreader.get_box_extent()
    origin = ncreader.get_box_origin()

    nparcels = ncreader.get_num_parcels(step)
    ttime = ncreader.get_dataset(step, "t")

    label = ""
    if coloring == "aspect-ratio":
        label = "aspect ratio"
        data = ncreader.get_aspect_ratio(step=step)
    elif coloring == "vol-distr":
        label = "volume distribution"
        data = ncreader.get_dataset(step=step, name="volume")
        data[data <= vmin] = 0.0
        data[data > vmin] = 1.0
    else:
        label = coloring
        data = ncreader.get_dataset(step=step, name=coloring)

    if title is None:
        title = (
            label
            + "\t\t\t\t time = %15.3f" % ttime
            + "\t\t\t\t #parcels = %10d" % nparcels
        )

    if no_title:
        title = None

    if graph is None:
        graph, _ = _get_bokeh_basic_graph(
            origin=origin, extent=extent, title=title, **kwargs
        )

    x, y, width, height, angle = ncreader.get_ellipses_for_bokeh(step)

    if norm:
        vmin = -1.0
        vmax = 1.0

        dmax = data.max()
        dmin = data.min()

        data[data > 0] = data[data > 0] / dmax
        data[data < 0] = data[data < 0] / (-dmin)


    source = ColumnDataSource(
        dict(x=x, y=y, width=width, height=height, angle=angle, fill_color=data)
    )

    font_size = bokeh_style["font.size"]
    text_font = bokeh_style["font.font"]

    mapper = None
    color_bar = None
    if coloring == "vol-distr":
        mapper = linear_cmap(
            field_name="fill_color", palette=["blue", "red"], low=0, high=1
        )
        # 5 August 2021
        # https://stackoverflow.com/questions/63344015/how-do-i-get-a-bokeh-colorbar-to-show-the-min-and-max-value
        ticker = FixedTicker(ticks=[0, 0.5, 1])
        color_bar = ColorBar(
            color_mapper=mapper["transform"],
            label_standoff=12,
            ticker=ticker,
            title_text_font_size=font_size,
            major_label_text_font=text_font,
            major_label_text_font_size=font_size,
        )
        # 5 August 2021
        # https://stackoverflow.com/questions/37173230/how-do-i-use-custom-labels-for-ticks-in-bokeh
        color_bar.major_label_overrides = {0: "0", 0.5: "Vmin", 1: "Vmax"}
    else:
        mapper = linear_cmap(
            field_name="fill_color", palette=palette, low=vmin, high=vmax
        )
        color_bar = ColorBar(
            color_mapper=mapper["transform"],
            label_standoff=12,
            title_text_font_size=font_size,
            major_label_text_font=text_font,
            major_label_text_font_size=font_size,
        )

    graph.ellipse(
        x="x",
        y="y",
        width="width",
        height="height",
        angle="angle",
        color=mapper,
        fill_alpha=fill_alpha,
        line_color=None,
        source=source,
    )

    if not no_colorbar:
        color_bar.width = color_bar_width
        graph.add_layout(color_bar, "right")

    return graph


def bokeh_plot(fname, step, show=False, fmt="png", coloring="vorticity", **kwargs):

    hybrid = kwargs.pop("hybrid", False)

    ncreader = nc_reader()

    ncreader.open(fname)

    if hybrid and not ncreader.is_field_file:
        if not ncreader.is_parcel_file:
            raise RuntimeError("Neither a field nor parcel file.")
        ncreader.close()
        fname = fname.replace("_" + str(step).zfill(10) + "_parcels.nc", "_fields.nc")
        ncreader.open(fname)

    nsteps = ncreader.get_num_steps()

    if step > nsteps - 1:
        raise ValueError("Step number exceeds limit of " + str(nsteps - 1) + ".")

    if step < 0:
        raise ValueError("Step number cannot be negative.")

    saveas = coloring + "_step_" + str(step).zfill(len(str(nsteps)))


    if ncreader.is_parcel_file:
        if coloring == "aspect-ratio":
            vmin = 1.0
            vmax = ncreader.get_global_attribute("lambda_max")
        elif coloring == "vol-distr":
            extent = ncreader.get_box_extent()
            ncells = ncreader.get_box_ncells()
            vcell = np.prod(extent / ncells)
            vmin = vcell / ncreader.get_global_attribute("min_vratio")
            vmax = vcell / ncreader.get_global_attribute("max_vratio")
        else:
            vmin, vmax = ncreader.get_dataset_min_max(coloring)

        saveas = "parcels_" + saveas
        graph = _bokeh_plot_parcels(ncreader, step, coloring, vmin, vmax, **kwargs)

    elif ncreader.is_field_file:
        saveas = "field_" + saveas
        col = coloring
        if col == "buoyancy":
            col = "total " + col
        vmin, vmax = ncreader.get_dataset_min_max(col)
        graph = _bokeh_plot_field(
            ncreader, step, col, vmin, vmax, hybrid=hybrid, **kwargs
        )

        if hybrid:
            fname = fname.replace("_fields.nc", "_" + str(step).zfill(10) + "_parcels.nc")
            ncreader.close()
            ncreader.open(fname)

            if not (nsteps == ncreader.get_num_steps()):
                raise RuntimeError(
                    "Field and parcel file do not have the same step count"
                )

            graph = _bokeh_plot_parcels(
                ncreader,
                step,
                coloring,
                vmin,
                vmax,
                graph=graph,
                no_colorbar=True,
                fill_alpha=1.0,
                **kwargs
            )
    else:
        raise RuntimeError("Neither a field nor parcel file.")

    ncreader.close()

    _bokeh_save(graph=graph, fname=saveas, fmt=fmt, show=show, **kwargs)
