import h5py
import os
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


    def get_mesh_extent()
        return np.array(self.h5file['mesh'].attrs['extent'])


    def get_mesh_grid()
        return np.array(self.h5file['mesh'].attrs['grid'])


    def get_mesh_origin()
        return np.array(self.h5file['mesh'].attrs['origin'])


    def get_parcel_dataset(step, name):
        s = 'step#' + str(step).zfill(10)
        if not name in self.h5file[s]['parcels'].keys():
            raise IOError("Parcel dataset '" + name "' unknown.")
        return np.array(self.h5file[s]['parcels'][name])


    def get_ellipses(self, step)
        position = self.get_parcel_dataset('position')
        B = self.get_parcel_dataset('B')
        V = self.get_parcel_dataset('volume')
        angle = get_parcel_dataset('orientation')

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
        B = self.get_parcel_dataset('B')
        V = self.get_parcel_dataset('volume')

        B22 = self._get_B22(B[0, :], B[1, :])
        return self._get_eigenvalue(B[0, :], B[1, :], B22)


    def _get_B22(self, B11, B12):
        return (1.0 + B12 ** 2) / B11


    def _get_eigenvalue(self, B11, B12, B22):
        return 0.5 * (B11 + B22) + np.sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)
