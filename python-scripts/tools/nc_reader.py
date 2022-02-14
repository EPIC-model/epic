import netCDF4 as nc
import os
import re
import numpy as np
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import EllipseCollection


class nc_reader:
    def __init__(self):
        self._ncfile = None

        # either a field or parcel file
        self._nctype = None

    def open(self, fname):
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")
        self._nctype = self.get_global_attribute("file_type")

        # if we read in a parcel file we pre-evaluate the number of
        # steps
        if self.is_parcel_file:
            self._loaded_step = -1
            basename = os.path.basename(fname)
            # 14 Feb 2022
            # https://stackoverflow.com/questions/15340582/python-extract-pattern-matches
            p = re.compile("(.*)_(\d*)_parcels.nc")
            result = p.search(basename)
            self._basename = result.group(1)
            self._dirname = os.path.dirname(fname)
            if self._dirname == '':
                self._dirname = '.'
            self._n_parcel_files = 0
            for ff in os.listdir(self._dirname):
                if self._basename in ff and '_parcels.nc' in ff:
                    self._n_parcel_files += 1



    def close(self):
        self._ncfile.close()

    @property
    def is_parcel_file(self):
        return self._nctype == "parcels"

    @property
    def is_parcel_stats_file(self):
        return self._nctype == "parcel_stats"

    @property
    def is_field_stats_file(self):
        return self._nctype == "field_stats"

    @property
    def is_field_file(self):
        return self._nctype == "fields"

    def get_num_steps(self):
        if self.is_parcel_file:
            return self._n_parcel_files
        else:
            return self._ncfile.dimensions['t'].size

    def get_box_extent(self):
        return self.get_global_attribute("extent")

    def get_box_ncells(self):
        return self.get_global_attribute("ncells")

    def get_box_origin(self):
        return self.get_global_attribute("origin")

    #TODO
    #def get_parcel_option(self, name):
        #if not name in self._ncfile["options"]["parcel"].attrs.keys():
            #raise IOError("Parcel options '" + name + "' unknown.")
        #opt = self._ncfile["options"]["parcel"].attrs[name][0]
        #if isinstance(opt, np.bytes_):
            ## 16 June 2021
            ## https://stackoverflow.com/questions/35576999/how-to-read-strings-from-hdf5-dataset-using-h5py
            #opt = opt.decode()
        #return opt


    def get_dataset(self, step, name, indices=None):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        if step > nsteps - 1:
            raise ValueError("Dataset has only steps 0 to " + str(nsteps - 1) + ".")

        if self.is_parcel_file:
            step = step + 1
            self._load_parcel_file(step)
            if indices is not None:
                return self._ncfile.variables[name][indices, ...]
            else:
                return np.array(self._ncfile.variables[name])
        else:
            if indices is not None:
                return self._ncfile.variables[name][step, ...][indices, ...]
            else:
                return np.array(self._ncfile.variables[name][step, ...])

    def get_dataset_min_max(self, name, indices=None):
        nsteps = self.get_num_steps()
        data = self.get_dataset(0, name, indices=indices)
        vmax = data.max()
        vmin = data.min()
        for step in range(1, nsteps):
            data = self.get_dataset(step, name, indices=indices)
            vmax = max(vmax, data.max())
            vmin = min(vmin, data.min())
        return vmin, vmax

    #TODO
    def get_global_attribute_names(self):
        return list(self._ncfile.ncattrs())

    #TODO
    def get_global_attribute(self, name):
        if not name in self._ncfile.ncattrs():
            raise IOError("Global attribute '" + name + "' unknown.")
        attr = self._ncfile.getncattr(name)
        if isinstance(attr, np.bytes_):
            attr = attr.decode()
        return attr

    #TODO
    def get_step_attribute_names(self):
        s = self._get_step_string(step)
        return list(self._ncfile[s].attrs.keys())

    #TODO
    def get_step_attribute(self, step, name):
        s = self._get_step_string(step)
        if not name in self._ncfile[s].attrs.keys():
            raise IOError("Step attribute '" + name + "' unknown.")
        val = self._ncfile[s].attrs[name]
        if isinstance(val, np.float64):
            return val
        return val[0]

    #TODO
    def get_diagnostic(self, name):
        s = self._get_step_string(0)
        if not name in self._ncfile[s].attrs.keys():
            raise IOError("Attribute '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        shape = np.array(self._ncfile[s].attrs[name].shape)
        shape[0] = nsteps
        data = np.zeros(shape)
        for step in range(nsteps):
            s = self._get_step_string(step)
            data[step] = np.array(self._ncfile[s].attrs[name])
        return data

    def get_num_parcels(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        self._load_parcel_file(step)
        return self._ncfile.dimensions['n_parcels'].size

    def get_ellipses(self, step, indices=None):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        position = self.get_dataset(step, "position", indices=indices)
        V = self.get_dataset(step, "volume", indices=indices)
        B = self.get_dataset(step, "B", indices=indices)

        B22 = self._get_B22(B[:, 0], B[:, 1], V)
        a2 = self._get_eigenvalue(B[:, 0], B[:, 1], B22)
        angle = self._get_angle(B[:, 0], B[:, 1], B22, a2)

        b2 = (V / np.pi) ** 2 / a2
        # 4 Feb 2022
        # https://matplotlib.org/stable/gallery/shapes_and_collections/ellipse_collection.html
        return EllipseCollection(widths=2 * np.sqrt(a2),
                                 heights=2 * np.sqrt(b2),
                                 angles=np.rad2deg(angle),
                                 units='xy',
                                 offsets=position)

    def get_ellipses_for_bokeh(self, step, indices=None):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        position = self.get_dataset(step, "position", indices=indices)
        V = self.get_dataset(step, "volume", indices=indices)
        B = self.get_dataset(step, "B", indices=indices)
        B22 = self._get_B22(B[:, 0], B[:, 1], V)
        a2 = self._get_eigenvalue(B[:, 0], B[:, 1], B22)
        angle = self._get_angle(B[:, 0], B[:, 1], B22, a2)
        b2 = (V / np.pi) ** 2 / a2
        return (
            position[:, 0],
            position[:, 1],
            2 * np.sqrt(a2[:]),
            2 * np.sqrt(b2[:]),
            angle[:],
        )

    def get_aspect_ratio(self, step, indices=None):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        V = self.get_dataset(step, "volume", indices=indices)
        B = self.get_dataset(step, "B", indices=indices)
        B22 = self._get_B22(B[:, 0], B[:, 1], V)
        a2 = self._get_eigenvalue(B[:, 0], B[:, 1], B22)
        return a2 / V * np.pi

    def _get_B22(self, B11, B12, V):
        return ((V / np.pi) ** 2 + B12 ** 2) / B11

    def _get_eigenvalue(self, B11, B12, B22):
        return 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)

    def _get_eigenvector(self, a2, B11, B12, B22):
        evec = np.array([a2 - B22, B12])

        for i in range(evec.shape[1]):
            if abs(evec[0, i]) + abs(evec[1, i]) == 0.0:
                if B11[i] > B22[i]:
                    evec[0, i] = evec[0, i] + np.finfo(np.float64).eps
                else:
                    evec[1, i] = evec[1, i] + np.finfo(np.float64).eps

        return evec / np.linalg.norm(evec, 2)

    def _get_angle(self, B11, B12, B22, a2=None):
        if a2 is None:
            a2 = self._get_eigenvalue(B11, B12, B22)
        evec = self._get_eigenvector(a2, B11, B12, B22)
        return np.arctan2(evec[1, :], evec[0, :])

    def _get_step_string(self, step):
        return str(step).zfill(10)

    def _load_parcel_file(self, step):
        if self._loaded_step == step:
            return
        self._loaded_step = step
        self._ncfile.close()
        s = self._get_step_string(step)
        fname = os.path.join(self._dirname, self._basename + '_' + s + '_parcels.nc')
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")
