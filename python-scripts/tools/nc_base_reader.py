import netCDF4 as nc
import numpy as np
from abc import ABC, abstractmethod
import os

class nc_base_reader(ABC):

    def __init__(self):
        self._ncfile = None

    def open(self, fname):
        if not os.path.exists(fname):
            raise IOError("File '" + fname + "' does not exist.")
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")

    def close(self):
        self._ncfile.close()

    def handle(self):
        return self._ncfile

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

    def has_global_attributes(self):
        return not self._ncfile.ncattrs() == []

    # 18 Feb 2022
    # https://stackoverflow.com/questions/8450472/how-to-print-a-string-at-a-fixed-width
    # 19 Feb 2022
    # https://stackoverflow.com/questions/873327/pythons-most-efficient-way-to-choose-longest-string-in-list
    def __str__(self):
        if self.has_global_attributes():
            print("=" * 80)
            # print global attributes
            print("GLOBAL ATTRIBUTES:")
            l = len(max(self._ncfile.ncattrs(), key=len))
            fmt = '{0: <' + str(l) + '}'
            for key in self._ncfile.ncattrs():
                print(fmt.format(key), "\t", self._ncfile.getncattr(key))
            print("-" * 80)

        print("DIMENSIONS:")

        for dim in self._ncfile.dimensions:
            print("    ", dim, "=", self._ncfile.dimensions[dim].size)

        print("-" * 80)

        print("VARIABLES:")
        # get first variable name
        name = list(self._ncfile.variables.keys())[0]

        if not self._ncfile.variables[name].ncattrs() == []:
            # get length of longest attribute string
            l = len(max(self._ncfile.variables[name].ncattrs(), key=len))
            fmt = '{0: <' + str(l) + '}'

        # print variables and their attributes
        for var in self._ncfile.variables:
            print("    ", var)
            for attr in self._ncfile.variables[var].ncattrs():
                print("\t", fmt.format(attr), "\t", self._ncfile.variables[var].getncattr(attr))
        print("=" * 80)
        return ""
