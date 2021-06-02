import numpy as np
from tools.h5_reader import H5Reader
import progressbar
from bokeh.io import curdoc, show, export_png
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, Grid, LinearAxis, Plot, ColorBar, LinearColorMapper
from bokeh.palettes import Viridis256
from bokeh.transform import linear_cmap
import matplotlib.pyplot as plt
from bokeh.io.export import get_screenshot_as_png
import os,subprocess,shutil

class BokehAnimation:
    def __init__(self):
        pass

    def create(self, fname, coloring='aspect-ratio', xmin=None, xmax=None, ymin=None, ymax=None):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        self.nsteps = self.h5reader.get_num_steps()
        self.extent = self.h5reader.get_mesh_extent()
        self.origin = self.h5reader.get_mesh_origin()

        os.mkdir("temp-movie")
        #file to save the model

        # instantiating the figure object
        left=self.origin[0]
        right=self.origin[0]+self.extent[0]
        bottom=self.origin[1]
        top=self.origin[1]+self.extent[1]
        if(xmin is not None):
           left=xmin
        if(xmax is not None):
           right=xmax
        if(ymin is not None):
           bottom=ymin
        if(ymax is not None):
           top=ymax
        pw=np.nanmin([1920,int(1080*(right-left)/(top-bottom))])
        ph=np.nanmin([1080,int(1920*(top-bottom)/(right-left))])
        output_file("temp-movie/movie-template.html")
        graph = figure(output_backend="webgl", title =coloring,plot_width=pw,plot_height=ph,match_aspect=True, x_range=(left, right), y_range=(bottom, top))
        step=0
        x,y,width,height,angle=self.h5reader.get_ellipses_for_bokeh(step)
        if (coloring=='aspect-ratio'):
            variable_of_interest=self.h5reader.get_aspect_ratio(step)
            variable_of_interest_min=1
            variable_of_interest_max=self.h5reader.get_parcel_info('lambda')[0]
        else:
            variable_of_interest=self.h5reader.get_parcel_dataset(step, coloring)
            variable_of_interest_min=np.nanmin(variable_of_interest)
            variable_of_interest_max=np.nanmax(variable_of_interest)
        source = ColumnDataSource(dict(x=x,y=y, width=width, height=height,angle=angle,fill_color=variable_of_interest))
        Viridis256_r = tuple(reversed(list(Viridis256)))
        mapper=linear_cmap(field_name='fill_color', palette=Viridis256_r, low=variable_of_interest_min, high=variable_of_interest_max)
        graph.ellipse(x='x', y='y', width='width', height='height',angle='angle',color=mapper,fill_alpha=0.75,line_color=None,source=source)
        color_bar = ColorBar(color_mapper=mapper['transform'], label_standoff=12)
        graph.add_layout(color_bar, 'right')
        for step in range(self.nsteps):
            x,y,width,height,angle=self.h5reader.get_ellipses_for_bokeh(step)
            if (coloring=='aspect-ratio'):
                variable_of_interest=self.h5reader.get_aspect_ratio(step)
            else:
                variable_of_interest=self.h5reader.get_parcel_dataset(step, coloring)
            new_source = dict(x=x,y=y, width=width, height=height,angle=angle,fill_color=variable_of_interest)
            source.data=new_source
            export_png(graph, filename = "temp-movie/movie.%05d.png"%step)

    def save(self, fname):
        ffmpeg_cmd1 = ['ffmpeg','-f','image2','-framerate','5','-i',os.path.join('temp-movie/movie.%05d.png'),'-c:','libx264','-pix_fmt','yuv444p','-vf','fps=5','-crf','20',fname]
        subprocess.call(ffmpeg_cmd1,shell=False)
        #shutil.rmtree("temp-movie")

