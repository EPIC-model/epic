import netCDF4 as nc
import numpy as np
from abc import ABC, abstractmethod

class nc_base_reader(ABC):

    def open(self, fname):
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")

    def close(self):
        self._ncfile.close()

    @abstractmethod
    def get_num_steps(self, **kwargs):
        pass

    @abstractmethod
    def get_box_extent(self, **kwargs):
        pass

    @abstractmethod
    def get_box_ncells(self, **kwargs):
        pass

    @abstractmethod
    def get_box_origin(self, **kwargs):
        pass

    @abstractmethod
    def get_dataset(self, step, name, **kwargs):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        if step > nsteps - 1:
            raise ValueError("Dataset has only steps 0 to " + str(nsteps - 1) + ".")

    @abstractmethod
    def get_axis(self, **kwargs):
        pass

    def get_all(self, name):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")
        return np.array(self._ncfile.variables[name])

    def get_meshgrid(self, **kwargs):
        x = self.get_axis('x', **kwargs)
        y = self.get_axis('y', **kwargs)
        z = self.get_axis('z', **kwargs)

        xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')
        assert np.all(xg[:, 0, 0] == x)
        assert np.all(yg[0, :, 0] == y)
        assert np.all(zg[0, 0, :] == z)

        return xg, yg, zg
