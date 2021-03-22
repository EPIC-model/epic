import h5py
import os
import numpy as np
from matplotlib.patches import Ellipse

class H5Reader:

    def __init__(self):
        self.h5file = None


    def open(self, fname):
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")
        self.h5file = h5py.File(fname, 'r')


    def close(self):
        self.h5file.close()

    def get_num_steps(self):
        return self.h5file['/'].attrs['nsteps'][0]

    def get_mesh_extent(self):
        return np.array(self.h5file['mesh'].attrs['extent'])


    def get_mesh_grid(self):
        return np.array(self.h5file['mesh'].attrs['grid'])


    def get_mesh_origin(self):
        return np.array(self.h5file['mesh'].attrs['origin'])

    def get_diagnostic(self, name):
        s = 'step#' + str(0).zfill(10)
        if not name in self.h5file[s]['diagnostics'].attrs.keys():
            raise IOError("Diagnostic '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        shape = np.array(self.h5file[s]['diagnostics'].attrs[name].shape)
        shape[0] = nsteps
        data = np.zeros(shape)
        for step in range(nsteps):
            s = 'step#' + str(step).zfill(10)
            data[step] = np.array(self.h5file[s]['diagnostics'].attrs[name])
        return data


    def get_parcel_dataset(self, step, name):
        s = 'step#' + str(step).zfill(10)
        if not name in self.h5file[s]['parcels'].keys():
            raise IOError("Parcel dataset '" + name + "' unknown.")
        return np.array(self.h5file[s]['parcels'][name])


    def get_step_attribute(self, step, name):
        s = 'step#' + str(step).zfill(10)
        if not name in self.h5file[s].attrs.keys():
            raise IOError("Step attribute '" + name + "' unknown.")
        return self.h5file[s].attrs[name]


    def get_ellipses(self, step):
        position = self.get_parcel_dataset(step, 'position')
        B = self.get_parcel_dataset(step, 'B')
        V = self.get_parcel_dataset(step, 'volume')
        angle = self.get_parcel_dataset(step, 'orientation')

        B22 = self._get_B22(B[0, :], B[1, :])
        lam = self._get_eigenvalue(B[0, :], B[1, :], B22)

        a2 = lam * V / np.pi
        b2 = V / np.pi / lam

        return [Ellipse(xy=position[:, i],
                        width=2 * np.sqrt(a2[i]),
                        height=2 * np.sqrt(b2[i]),
                        angle=np.rad2deg(angle[i]))
                for i in range(B.shape[1])]


    def get_aspect_ratio(self, step):
        B = self.get_parcel_dataset(step, 'B')
        V = self.get_parcel_dataset(step, 'volume')

        B22 = self._get_B22(B[0, :], B[1, :])
        return self._get_eigenvalue(B[0, :], B[1, :], B22)


    def _get_B22(self, B11, B12):
        return (1.0 + B12 ** 2) / B11


    def _get_eigenvalue(self, B11, B12, B22):
        return 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)
