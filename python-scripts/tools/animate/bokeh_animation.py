import numpy as np
from tools.nc_reader import nc_reader
import progressbar
from bokeh.io import export_png
import matplotlib.pyplot as plt
import os, subprocess, shutil
from tools.bokeh_plots import _bokeh_plot_parcels


class BokehAnimation:
    def __init__(self):
        pass

    def create(self, fname, coloring="aspect-ratio", **kwargs):
        self.ncreader = nc_reader()

        self.ncreader.open(fname)

        if not self.ncreader.is_parcel_file:
            raise IOError("Not a parcel output file.")

        nsteps = self.ncreader.get_num_steps()
        os.mkdir("temp-movie")

        tmin = kwargs.pop("tmin", None)
        tmax = kwargs.pop("tmax", None)

        if tmin is None:
            tmin = 0.0

        if tmax is None:
            t1 = self.ncreader.get_dataset(nsteps - 2, "t")
            t2 = self.ncreader.get_dataset(nsteps - 3, "t")
            dt = t1 - t2
            tmax = t1 + 10 * dt

        bar = progressbar.ProgressBar(maxval=nsteps).start()
        i = 0
        for step in range(0, nsteps):

            t = self.ncreader.get_dataset(step, "t")
            if t < tmin:
                continue

            if t > tmax:
                break

            if coloring == "aspect-ratio":
                vmin = 1.0
                vmax = self.ncreader.get_global_attribute("lambda_max")
            else:
                vmin, vmax = self.ncreader.get_dataset_min_max(coloring)

            graph = _bokeh_plot_parcels(
                self.ncreader, step, coloring, vmin, vmax, **kwargs
            )

            export_png(graph, filename="temp-movie/movie.%05d.png" % i)
            bar.update(step + 1)
            i = i + 1
        bar.finish()

    def save(self, fname):
        ffmpeg_cmd1 = [
            "ffmpeg",
            "-f",
            "image2",
            "-framerate",
            "5",
            "-i",
            os.path.join("temp-movie/movie.%05d.png"),
            "-c:",
            "libx264",
            "-pix_fmt",
            "yuv444p",
            "-vf",
            "fps=5",
            "-crf",
            "20",
            fname,
        ]
        subprocess.call(ffmpeg_cmd1, shell=False)
        shutil.rmtree("temp-movie")
