import netCDF4 as nc
import os
import numpy as np


class nc_parcels:
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
        self._nparcels = 0


    def add_dataset(self, name, values, dtype='f8', **kwargs):
        """
        Add a parcel attribute dataset.

        Parameters
        ----------
        name: str
            The field name.
        values: np.ndarray
            The field data.
        """

        values = np.asarray(values)

        shape = np.shape(values)
        if len(shape) > 1:
            RuntimeError("Shape must be of 1-dimensional.")

        if self._nparcels == 0:
            self._nparcels = len(values)
            # add dimension
            self._ncfile.createDimension(dimname="n_parcels", size=self._nparcels)
            self._ncfile.setncattr('file_type', 'parcels')
            self._ncfile.setncattr('t', 0.0)

        var = self._ncfile.createVariable(varname=name,
                                          datatype=dtype,
                                          dimensions=("n_parcels"))
        var[:] = values[:]

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
