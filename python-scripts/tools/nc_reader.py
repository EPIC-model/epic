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

        self._is_compressible = False

        self._derived_fields = [
            'vorticity_magnitude',
            'helicity',
            'enstrophy',
            'cross_helicity_magnitude',
            'kinetic_energy'
        ]

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

    def get_axis(self, name):
        axis = self.get_all(name)
        if name == 'x' or name == 'y':
            # copy periodic grid point
            axis = np.append(axis, abs(axis[0]))
        return

    def get_meshgrid(self):
        x = self.get_axis('x')
        y = self.get_axis('y')
        z = self.get_axis('z')

        xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')

        # 13 July 2022
        # https://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d
        assert np.all(xg[:, 0, 0] == x)
        assert np.all(yg[0, :, 0] == y)
        assert np.all(zg[0, 0, :] == z)

        return xg, yg, zg

    def get_all(self, name):
        if self.is_parcel_file:
            raise IOError("This function is not availble for parcel files.")

        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        return np.array(self._ncfile.variables[name])

    def get_dataset(self, step, name, indices=None, copy_periodic=True):

        if name in self._derived_fields:
            return self._get_derived_dataset(step, name, copy_periodic)

        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        nsteps = self.get_num_steps()
        if step > nsteps - 1:
            raise ValueError("Dataset has only steps 0 to " + str(nsteps - 1) + ".")

        if self.is_parcel_file and name == 't':
            # parcel files store the time as a global attribute
            self._load_parcel_file(step)
            return np.array(self._ncfile.variables[name]).squeeze()

        if self.is_parcel_file:
            self._load_parcel_file(step)
            if indices is not None:
                return np.array(self._ncfile.variables[name]).squeeze()[indices, ...]
            else:
                return np.array(self._ncfile.variables[name]).squeeze()
        else:
            if indices is not None:
                return np.array(self._ncfile.variables[name][step, ...]).squeeze()[indices, ...]
            else:
                fdata = np.array(self._ncfile.variables[name][step, ...]).squeeze()

                if copy_periodic:
                    fdata = self._copy_periodic_layers(fdata)

                # change ordering from (z, y, x) to (x, y, z)
                fdata = np.transpose(fdata, axes=[2, 1, 0])

                return fdata

    def _get_derived_dataset(self, step, name, copy_periodic):
        if name == 'vorticity_magnitude':
            x_vor = self.get_dataset(step=step, name='x_vorticity', copy_periodic=copy_periodic)
            y_vor = self.get_dataset(step=step, name='y_vorticity', copy_periodic=copy_periodic)
            z_vor = self.get_dataset(step=step, name='z_vorticity', copy_periodic=copy_periodic)
            return np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)
        if name == 'helicity':
            u = self.get_dataset(step=step, name='x_velocity', copy_periodic=copy_periodic)
            v = self.get_dataset(step=step, name='y_velocity', copy_periodic=copy_periodic)
            w = self.get_dataset(step=step, name='z_velocity', copy_periodic=copy_periodic)
            xi = self.get_dataset(step=step, name='x_vorticity', copy_periodic=copy_periodic)
            eta = self.get_dataset(step=step, name='y_vorticity', copy_periodic=copy_periodic)
            zeta = self.get_dataset(step=step, name='z_vorticity', copy_periodic=copy_periodic)
            return u * xi + v * eta + w * zeta
        if name == 'cross_helicity_magnitude':
            u = self.get_dataset(step=step, name='x_velocity', copy_periodic=copy_periodic)
            nx, ny, nz = u.shape
            uvec = np.zeros((nx, ny, nz, 3))
            ovec = np.zeros((nx, ny, nz, 3))
            uvec[:, :, :, 0] = u
            uvec[:, :, :, 1] = self.get_dataset(step=step, name='y_velocity',
                                                copy_periodic=copy_periodic)
            uvec[:, :, :, 2] = self.get_dataset(step=step, name='z_velocity',
                                                copy_periodic=copy_periodic)
            ovec[:, :, :, 0] = self.get_dataset(step=step, name='x_vorticity',
                                                copy_periodic=copy_periodic)
            ovec[:, :, :, 1] = self.get_dataset(step=step, name='y_vorticity',
                                                copy_periodic=copy_periodic)
            ovec[:, :, :, 2] = self.get_dataset(step=step, name='z_vorticity',
                                                copy_periodic=copy_periodic)
            ch = np.cross(uvec, ovec)
            x_ch = ch[:, :, :, 0]
            y_ch = ch[:, :, :, 1]
            z_ch = ch[:, :, :, 2]
            return np.sqrt(x_ch ** 2 + y_ch  ** 2 + z_ch ** 2)

        if name == 'kinetic_energy':
            u = self.get_dataset(step=step, name='x_velocity', copy_periodic=copy_periodic)
            v = self.get_dataset(step=step, name='y_velocity', copy_periodic=copy_periodic)
            w = self.get_dataset(step=step, name='z_velocity', copy_periodic=copy_periodic)
            return 0.5 * (u ** 2 + v ** 2 + w ** 2)
        if name == 'enstrophy':
            xi = self.get_dataset(step=step, name='x_vorticity', copy_periodic=copy_periodic)
            eta = self.get_dataset(step=step, name='y_vorticity', copy_periodic=copy_periodic)
            zeta = self.get_dataset(step=step, name='z_vorticity', copy_periodic=copy_periodic)
            return 0.5 * (xi ** 2 + eta ** 2 + zeta ** 2)

    def get_dataset_attribute(self, name, attr):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")

        if not attr in self._ncfile.variables[name].ncattrs():
            raise IOError("Dataset attribute '" + name + "' unknown.")

        return self._ncfile.variables[name].getncattr(attr)

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

    def get_global_attribute_names(self):
        return list(self._ncfile.ncattrs())

    def get_global_attribute(self, name):
        if not name in self._ncfile.ncattrs():
            raise IOError("Global attribute '" + name + "' unknown.")
        attr = self._ncfile.getncattr(name)
        if isinstance(attr, np.bytes_):
            attr = attr.decode()
        return attr

    def get_diagnostic(self, name):
        if not name in self._ncfile.variables.keys():
            raise IOError("Dataset '" + name + "' unknown.")
        return np.array(self._ncfile.variables[name])

    def get_num_parcels(self, step):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        self._load_parcel_file(step)
        return self._ncfile.dimensions['n_parcels'].size

    def get_ellipses(self, step, indices=None):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        x_pos = self.get_dataset(step, "x_position", indices=indices)
        z_pos = self.get_dataset(step, "z_position", indices=indices)
        position = np.empty((len(x_pos), 2))
        position[:, 0] = x_pos
        position[:, 1] = z_pos
        V = self.get_dataset(step, "volume", indices=indices)
        B11 = self.get_dataset(step, "B11", indices=indices)
        B12 = self.get_dataset(step, "B12", indices=indices)
        if self._is_compressible:
            B22 = self.get_dataset(step, "B22", indices=indices)
        else:
            B22 = self._get_B22(B11, B12, V)
        a2 = self._get_eigenvalue(B11, B12, B22)
        angle = self._get_angle(B11, B12, B22, a2)

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
        x_pos = self.get_dataset(step, "x_position", indices=indices)
        z_pos = self.get_dataset(step, "z_position", indices=indices)
        V = self.get_dataset(step, "volume", indices=indices)
        B11 = self.get_dataset(step, "B11", indices=indices)
        B12 = self.get_dataset(step, "B12", indices=indices)
        if self._is_compressible:
            B22 = self.get_dataset(step, "B22", indices=indices)
        else:
            B22 = self._get_B22(B11, B12, V)
        a2 = self._get_eigenvalue(B11, B12, B22)
        angle = self._get_angle(B11, B12, B22, a2)
        #print("bye")
        #exit()
        b2 = (V / np.pi) ** 2 / a2
        #print(x_pos.shape, z_pos.shape, a2.shape, b2.shape, angle.shape)
        #exit()
        return (
            x_pos[:],
            z_pos[:],
            2 * np.sqrt(a2[:]),
            2 * np.sqrt(b2[:]),
            angle[:],
        )

    def get_aspect_ratio(self, step, indices=None):
        if not self.is_parcel_file:
            raise IOError("Not a parcel output file.")
        V = self.get_dataset(step, "volume", indices=indices)
        B11 = self.get_dataset(step, "B11", indices=indices)
        B12 = self.get_dataset(step, "B12", indices=indices)
        if self._is_compressible:
            B22 = self.get_dataset(step, "B22", indices=indices)
        else:
            B22 = self._get_B22(B11, B12, V)
        a2 = self._get_eigenvalue(B11, B12, B22)
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
        step = step + 1
        if self._loaded_step == step:
            return
        self._loaded_step = step
        self._ncfile.close()
        s = self._get_step_string(step)
        fname = os.path.join(self._dirname, self._basename + '_' + s + '_parcels.nc')
        self._ncfile = nc.Dataset(fname, "r", format="NETCDF4")

    # 18 Feb 2022
    # https://stackoverflow.com/questions/8450472/how-to-print-a-string-at-a-fixed-width
    # 19 Feb 2022
    # https://stackoverflow.com/questions/873327/pythons-most-efficient-way-to-choose-longest-string-in-list
    def __str__(self):
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

    def _copy_periodic_layers(self, field):
        nz, ny, nx = field.shape
        field_copy = np.empty((nz, ny+1, nx+1))
        field_copy[:, 0:ny, 0:nx] = field.copy()
        field_copy[:, ny, :] = field_copy[:, 0, :]
        field_copy[:, :, nx] = field_copy[:, :, 0]
        return field_copy
