import h5py
import os
import numpy as np
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import EllipseCollection


class H5Reader:
    def __init__(self):
        self._h5file = None

        # either a field or parcel file
        self._h5type = None

    def open(self, fname):
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")
        self._h5file = h5py.File(fname, "r")
        self._h5type = self.get_global_attribute("output_type")

    def close(self):
        self._h5file.close()

    @property
    def is_parcel_file(self):
        return self._h5type == "parcels"

    @property
    def is_parcel_stats_file(self):
        return self._h5type == "parcel diagnostics"

    @property
    def is_field_stats_file(self):
        return self._h5type == "field diagnostics"

    @property
    def is_field_file(self):
        return self._h5type == "fields"

    def get_num_steps(self):
        return self._h5file.attrs["nsteps"][0]

    def get_box_extent(self):
        return np.array(self._h5file["box"].attrs["extent"])

    def get_box_ncells(self):
        return np.array(self._h5file["box"].attrs["ncells"])

    def get_parcel_option(self, name):
        if not name in self._h5file["options"]["parcel"].attrs.keys():
            raise IOError("Parcel options '" + name + "' unknown.")
        opt = self._h5file["options"]["parcel"].attrs[name][0]
        if isinstance(opt, np.bytes_):
            # 16 June 2021
            # https://stackoverflow.com/questions/35576999/how-to-read-strings-from-hdf5-dataset-using-h5py
            opt = opt.decode()
        return opt

    def get_box_origin(self):
        return np.array(self._h5file["box"].attrs["origin"])

    def get_dataset(self, step, name, indices=None):
        s = self._get_step_string(step)
        if not name in self._h5file[s].keys():
            raise IOError("Dataset '" + name + "' unknown.")
        if indices is not None:
            return np.array(self._h5file[s][name])[indices, ...]
        else:
            return np.array(self._h5file[s][name])

    def get_dataset_min_max(self, name):
        nsteps = self.get_num_steps()
        data = self.get_dataset(0, name)
        vmax = data.max()
        vmin = data.min()
        for step in range(1, nsteps):
            data = self.get_dataset(step, name)
            vmax = max(vmax, data.max())
            vmin = min(vmin, data.min())
        return vmin, vmax

    def get_global_attribute_names(self):
        return list(self._h5file.attrs.keys())

    def get_global_attribute(self, name):
        if not name in self._h5file.attrs.keys():
            raise IOError("Global attribute '" + name + "' unknown.")
        attr = self._h5file.attrs[name][0]
        if isinstance(attr, np.bytes_):
            attr = attr.decode()
        return attr

    def get_step_attribute_names(self):
        s = self._get_step_string(step)
        return list(self._h5file[s].attrs.keys())

    def get_step_attribute(self, step, name):
        s = self._get_step_string(step)
        if not name in self._h5file[s].attrs.keys():
            raise IOError("Step attribute '" + name + "' unknown.")
        val = self._h5file[s].attrs[name]
        if isinstance(val, np.float64):
            return val
        return val[0]

    def get_diagnostic_names(self):
        return self.get_step_attribute_names()

    def get_diagnostic(self, name):
        s = self._get_step_string(0)
        if not name in self._h5file[s].attrs.keys():
            raise IOError("Attribute '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        shape = np.array(self._h5file[s].attrs[name].shape)
        shape[0] = nsteps
        data = np.zeros(shape)
        for step in range(nsteps):
            s = self._get_step_string(step)
            data[step] = np.array(self._h5file[s].attrs[name])
        return data

    def get_num_parcels(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        return self.get_step_attribute(step, "num parcel")

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
        #return [
            #Ellipse(
                #xy=position[i, :],
                #width=2 * np.sqrt(a2[i]),
                #height=2 * np.sqrt(b2[i]),
                #angle=np.rad2deg(angle[i]),
            #)
            #for i in range(len(V))
        #]

    def get_ellipses_for_bokeh(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        position = self.get_dataset(step, "position")
        V = self.get_dataset(step, "volume")
        B = self.get_dataset(step, "B")
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

    def get_aspect_ratio(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        V = self.get_dataset(step, "volume")
        B = self.get_dataset(step, "B")
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
        return "step#" + str(step).zfill(10)
