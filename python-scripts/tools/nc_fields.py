import netCDF4 as nc
import os
import numpy as np
from tools.nc_utils import write_nc_info, write_nc_parameters


class nc_fields:
    def __init__(self):
        self._ncfile = None

    def open(self, fname):
        """
        Open an NetCDF file.

        Parameters
        ----------
        fname : str
            Filename with extension.
        """
        self._ncfile = nc.Dataset(fname, "w", format="NETCDF4")

        write_nc_info(ncfile=self._ncfile, file_type='fields')

        self._ndims = 0

        self._physical_quantities = {}

        self._parameters = {}

        self._time = 0.0

    def set_time(self, time):
        self._time = time

    def add_physical_quantity(self, key, value):
        self._physical_quantities[key] = value

    def add_parameter(self, key, value):
        self._parameters[key] = value

    def add_axis(self, axis, values):
        if axis == 'x' or axis == 'y' or axis == 'z':
            var = self._ncfile.createVariable(varname=axis,
                                              datatype='f8',                                            
                                              dimensions=(axis))
            var[:] = values[:]

    def add_field(self, name, values, dtype='f8', **kwargs):
        """
        Add a field dataset.

        Parameters
        ----------
        name: str
            The field name.
        values: np.ndarray
            The field data.
        """

        values = np.asarray(values)

        if self._ndims == 0:
            shape = np.shape(values)

            # add dimensions
            self._ncfile.createDimension(dimname="t", size=None)
            if len(shape) == 2:
                self._ncfile.createDimension(dimname="z", size=shape[0])
                self._ncfile.createDimension(dimname="x", size=shape[1])
                self._ndims = 2
            elif len(shape) == 3:
                self._ncfile.createDimension(dimname="z", size=shape[0])
                self._ncfile.createDimension(dimname="y", size=shape[1])
                self._ncfile.createDimension(dimname="x", size=shape[2])
                self._ndims = 3
            else:
                RuntimeError("Shape must be of 2 or 3 dimensions")

            time = self._ncfile.createVariable(varname='t',
                                               datatype='f8',
                                               dimensions=('t'))
            time[0] = self._time


        if self._ndims == 2:
            var = self._ncfile.createVariable(varname=name,
                                              datatype=dtype,
                                              dimensions=('t', 'z', 'x'))
            var[0, :, :] = values[:, :]
        else:
            var = self._ncfile.createVariable(varname=name,
                                              datatype=dtype,
                                              dimensions=('t', 'z', 'y', 'x'))
            var[0, :, :, :] = values[:, :, :]

        unit = kwargs.pop('unit', '')
        if unit:
            var.units = unit

        standard_name = kwargs.pop('standard_name', '')
        if standard_name:
            var.standard_name = standard_name

        long_name = kwargs.pop('long_name', '')
        if long_name:
            var.long_name = long_name

    def close(self):
        if not self._physical_quantities == {}:
            write_nc_parameters(self._ncfile, 'physical_quantities',
                                self._physical_quantities)

        if not self._parameters == {}:
            write_nc_parameters(self._ncfile, 'parameters',
                                self._parameters)

        self._ncfile.close()


    def add_box(self, origin, extent, ncells):
        """
        Box dictionary:

        Parameters
        ----------
        origin : np.array of floats (length 2 or 3)
            The origin of the domain (x, y, z).
        extent : np.array of floats (length 2 or 3)
            The extent of the box (x, y, z).
        ncells : np.array of floats (length 2 or 3)
            The number of cells per dimension (x, y, z).
        """
        origin = np.asarray(origin, dtype=np.float64)
        extent = np.asarray(extent, dtype=np.float64)
        ncells = np.asarray(ncells, dtype=np.int32)

        l = len(origin)

        if l < 2 or l > 3:
            raise RuntimeError("Array 'origin' must have length 2 or 3.")

        if not len(extent) == l:
            raise RuntimeError("Array 'extent' must have length 2 or 3.")

        if not len(ncells) == l:
            raise RuntimeError("Array 'ncells' must have length 2 or 3.")

        self._ncfile.setncattr(name="origin", value=origin)
        self._ncfile.setncattr(name="extent", value=extent)
        self._ncfile.setncattr(name="ncells", value=ncells)
