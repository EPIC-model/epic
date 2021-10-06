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
        self._h5file = h5py.File(fname, "w")

    def add_fields(self, origin, extent, fields):
        """
        Write all fields provided by the dictionary.

        Parameters
        ----------
        origin : tuple of floats (length 2)
            The origin of the domain (horizontal, vertical).
        extent : tuple of floats (length 2)
            The extent of the box (horizontal, vertical).
        fields : dict
            All fields to be added. They keys are the field
            names and the values are the data (horizontal, vertical).
        """
        if not isinstance(fields, dict):
            TypeError("Argument not a dictionary.")

        shape = None

        for key, values in fields.items():
            if shape is None:
                shape = np.shape(values)
            elif not shape == np.shape(values):
                raise RuntimeError("Field '" + key + "' must have shape " + str(shape))

            self._add_dataset(name=key, data=values)

        ncells = (shape[0], shape[1] - 1)
        self._add_box(origin, extent, ncells)

    def close(self):
        self._h5file.close()

    def _add_box(self, origin, extent, ncells):
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

        group = self._h5file.create_group("box")
        group.attrs.create(name="origin", data=origin, dtype="f8")
        group.attrs.create(name="extent", data=extent, dtype="f8")
        group.attrs.create(name="ncells", data=ncells, dtype="i4")

    def _add_dataset(self, name, data, dtype="f8"):
        """
        Add a field dataset. The horizontal dimension is the
        leading dimension (horizontal x vertical).
        """
        data = np.asarray(data)
        self._h5file.create_dataset(name=name, dtype=dtype, data=data)
