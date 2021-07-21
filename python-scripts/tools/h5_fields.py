import h5py
import os
import numpy as np

class h5fields:

    def __init__(self):
        self._h5file = None


    def open(self, fname):
        """
        Open an HDF5 file.

        Parameters
        ----------
        fname : str
            Filename with extension.
        """
        self._h5file = h5py.File(fname, 'w')


    def add_box(self, origin, extent, ncells):
        """
        Box dictionary:

        Parameters
        ----------
        origin : tuple of floats (length 2)
            The origin of the domain (horizontal, vertical).
        extent : tuple of floats (length 2)
            The extent of the box (horizontal, vertical).
        ncells : tuple of ints (length 2)
            The number of cells per dimension (horizontal, vertical).
        """
        origin = tuple(origin)
        extent = tuple(extent)
        ncells = tuple(ncells)

        if not len(origin) == 2:
            raise RuntimeError("Tuple 'origin' must have length 2.")

        if not len(extent) == 2:
            raise RuntimeError("Tuple 'extent' must have length 2.")

        if not len(ncells) == 2:
            raise RuntimeError("Tuple 'ncells' must have length 2.")

        group = self._h5file.create_group('box')
        group.attrs['origin'] = origin
        group.attrs['extent'] = extent
        group.attrs['ncells'] = ncells


    def add_dataset(self, name, data, dtype='f'):
        """
        Add a field dataset. The horizontal dimension is the
        leading dimension (horizontal x vertical).
        """
        data = np.asarray(data)
        data = data.transpose()
        self._h5file.create_dataset(name=name, dtype=dtype, data=data)


    def close(self):
        self._h5file.close()
