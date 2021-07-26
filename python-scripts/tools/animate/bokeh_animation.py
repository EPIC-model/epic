import numpy as np
from tools.h5_reader import H5Reader
import progressbar
from bokeh.io import export_png
import matplotlib.pyplot as plt
import os, subprocess, shutil
from tools.bokeh_plots import _bokeh_plot_parcels

class BokehAnimation:
    def __init__(self):
        pass

    def create(self, fname, coloring='aspect-ratio', **kwargs):
        self.h5reader = H5Reader()

        self.h5reader.open(fname)

        if not self.h5reader.is_parcel_file:
            raise IOError('Not a parcel output file.')

        nsteps = self.h5reader.get_num_steps()
        os.mkdir("temp-movie")

        bar = progressbar.ProgressBar(maxval=nsteps).start()
        for step in range(nsteps):
            if coloring == 'aspect-ratio':
                vmin = 1.0
                vmax = self.h5reader.get_parcel_option('lambda')
            else:
                vmin, vmax = self.h5reader.get_dataset_min_max(coloring)

            graph = _bokeh_plot_parcels(self.h5reader, step, coloring, vmin, vmax, **kwargs)

            export_png(graph, filename = "temp-movie/movie.%05d.png"%step)
            bar.update(step+1)
        bar.finish()

    def save(self, fname):
        ffmpeg_cmd1 = ['ffmpeg','-f','image2','-framerate','5','-i',
        os.path.join('temp-movie/movie.%05d.png'),'-c:','libx264',
        '-pix_fmt','yuv444p','-vf','fps=5','-crf','20',fname]
        subprocess.call(ffmpeg_cmd1,shell=False)
        shutil.rmtree("temp-movie")

