import numpy as np
from tools.h5_reader import H5Reader
import progressbar
from bokeh.io import curdoc, show, export_png
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, Grid, LinearAxis, Plot, ColorBar, LinearColorMapper
from bokeh.palettes import Viridis
from bokeh.transform import linear_cmap
import matplotlib.pyplot as plt
from bokeh.io.export import get_screenshot_as_png
import os,subprocess,shutil

class BokehAnimation:
    def __init__(self):
        pass

    def create(self, fname, coloring='aspect-ratio'):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        self.nsteps = self.h5reader.get_num_steps()
        self.extent = self.h5reader.get_mesh_extent()
        self.origin = self.h5reader.get_mesh_origin()

        os.mkdir("temp-movie")
        #file to save the model
        pw=np.nanmin([1920,int(1080*self.extent[0]/self.extent[1])])
        ph=np.nanmin([1080,int(1920*self.extent[1]/self.extent[0])])
        # instantiating the figure object
        left=self.origin[0]
        right=self.origin[0]+self.extent[0]
        bottom=self.origin[1]
        top=self.origin[1]+self.extent[1]
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
        variable_scale=(variable_of_interest-variable_of_interest_min)/(variable_of_interest_max-variable_of_interest_min)
        cc_map=plt.get_cmap('viridis')
        fill_color=["#{:02x}{:02x}{:02x}".format(int(i[0]*255),int(i[1]*255),int(i[2]*255)) for i in cc_map(variable_scale)]
        source = ColumnDataSource(dict(x=x,y=y, width=width, height=height,angle=angle,fill_color=fill_color))
        graph.ellipse(x='x', y='y', width='width', height='height',angle='angle',fill_color='fill_color',fill_alpha=0.5,line_color=None,source=source)
        color_mapper = LinearColorMapper(palette="Viridis256", low=variable_of_interest_min, high=variable_of_interest_max)
        color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12)
        graph.add_layout(color_bar, 'right')
        for step in range(self.nsteps):
            x,y,width,height,angle=self.h5reader.get_ellipses_for_bokeh(step)
            if (coloring=='aspect-ratio'):
                variable_of_interest=self.h5reader.get_aspect_ratio(step)
            else:
                variable_of_interest=self.h5reader.get_parcel_dataset(step, coloring)
            variable_scale=(variable_of_interest-variable_of_interest_min)/(variable_of_interest_max-variable_of_interest_min)
            fill_color=["#{:02x}{:02x}{:02x}".format(int(i[0]*255),int(i[1]*255),int(i[2]*255)) for i in cc_map(variable_scale)]
            new_source = dict(x=x,y=y, width=width, height=height,angle=angle,fill_color=fill_color)
            source.data=new_source
            export_png(graph, filename = "temp-movie/movie.%05d.png"%step)

    def save(self, fname):
        ffmpeg_cmd1 = ['ffmpeg','-f','image2','-framerate','10','-i',os.path.join('temp-movie/movie.%05d.png'),'-c:v','mpeg4','-mbd','2','-trellis','2','-qscale:v','1',
'-pre_dia_size','4','-dia_size','4','-precmp','4','-cmp','4','-subcmp','4','-preme','2','-vf','fps=25','-b','3000k','-pass','1','-f','rawvideo','-y','/dev/null']
        ffmpeg_cmd2 = ['ffmpeg','-y','-f','image2','-framerate','10','-i',os.path.join('temp-movie/movie.%05d.png'),'-c:v','mpeg4','-mbd','2','-trellis','2','-qscale:v','1',
'-pre_dia_size','4','-dia_size','4','-precmp','4','-cmp','4','-subcmp','4','-preme','2','-vf','fps=25','-b','3000k','-pass','2', fname]
        subprocess.call(ffmpeg_cmd1,shell=False)
        subprocess.call(ffmpeg_cmd2,shell=False)
        os.remove("ffmpeg2pass-0.log")
        shutil.rmtree("temp-movie")

