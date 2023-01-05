from nc_base_reader import nc_base_reader
import re
import numpy as np
from matplotlib.patches import Ellipse, Circle
from matplotlib.collections import EllipseCollection
from geometry import ellipsoid, ellipse, xy_plane, xz_plane, yz_plane


class nc_reader(nc_base_reader):
    def __init__(self):
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
        super().open(fname)

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

    @property
    def is_three_dimensional(self):
        return len(self.get_box_ncells()) == 3

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
        if not name in ['x', 'y', 'z']:
            raise ValueError("No axis called '" + name + "'.")
        axis = self.get_all(name)
        if name == 'x' or name == 'y':
            axis = np.append(axis, abs(axis[0]))
        return axis

    def get_all(self, name):
        if self.is_parcel_file:
            raise IOError("This function is not availble for parcel files.")
        super().get_all(name)

    def get_dataset(self, step, name, indices=None, copy_periodic=True):
        super().get_dataset(step, name)

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

                if self.is_three_dimensional:
                    # change ordering from (z, y, x) to (x, y, z)
                    fdata = np.transpose(fdata, axes=[2, 1, 0])
                else:
                    fdata = np.transpose(fdata, axes=[1, 0])

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

    def get_physical_quantity(self, name):
        if not name in self._ncfile['physical_quantities'].ncattrs():
            raise IOError("Physical quantity '" + name + "' unknown.")
        attr = self._ncfile['physical_quantities'].getncattr(name)
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
        if self.is_three_dimensional:
            raise IOError("Not a 2-dimensional dataset.")

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
        if self.is_three_dimensional:
            raise IOError("Not a 2-dimensional dataset.")

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
        if self.is_three_dimensional:
            raise IOError("Not a 2-dimensional dataset.")

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

    def _get_B33(self, B11, B12, B13, B22, B23, V):
        return ((0.75 * V / np.pi) ** 2 - B13 * (B12 * B23 - B13 * B22) \
                                        + B11 * B23 ** 2                \
                                        - B12 * B13 * B23)              \
                / (B11 * B22 - B12 ** 2)

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

    def _copy_periodic_layers(self, field):
        if self.is_three_dimensional:
            nz, ny, nx = field.shape
            field_copy = np.empty((nz, ny+1, nx+1))
            field_copy[:, 0:ny, 0:nx] = field.copy()
            field_copy[:, ny, :] = field_copy[:, 0, :]
            field_copy[:, :, nx] = field_copy[:, :, 0]
        else:
            nz, nx = field.shape
            field_copy = np.empty((nz, nx+1))
            field_copy[:, 0:nx] = field.copy()
            field_copy[:, nx] = field_copy[:, 0]
        return field_copy

    def get_intersection_ellipses(self, step, plane, loc):
        """
        Calculates the ellipses from all ellipsoids intersecting
        with the provided xy-, xz- or yz-plane.

        Parameters
        ----------
        plane   'xy', 'xz' or 'yz'
        loc     grid point
        """
        if plane not in ['xy', 'xz', 'yz']:
            raise ValueError("Specified plane '", plane "' not available.")

        origin = self.get_origin()
        extent = self.get_box_extent()
        ncells = self.get_box_ncells()

        dx = extent / ncells

        var = {'yz': 0, 'xz': 1, 'xy': 2}
        dim = {0: 'x', 1: 'y', 2: 'z'}

        # calculate the position of the plane
        j = var[plane]
        p = origin[j] + dx[j] * loc

        if plane == 'xy':
            pl = xy_plane(z=p)
        elif plane ==  'xz':
            pl = xz_plane(y=p)
        else:
            pl = yz_plane(x=p)

        # lower and upper position bounds
        lo = origin[j] + dx[j] * (loc - 1)
        hi = origin[j] + dx[j] * (loc + 1)

        # get indices of parcels satisfying lo < pos < hi
        pos = self.get_dataset(step, name=dim[j] + '_position')
        indices = np.where((pos > lo) & (pos < hi))[0]
        pos = None

        B11 = self.get_dataset(step, name='B11',        indices=indices)
        B12 = self.get_dataset(step, name='B12',        indices=indices)
        B13 = self.get_dataset(step, name='B13',        indices=indices)
        B22 = self.get_dataset(step, name='B22',        indices=indices)
        B23 = self.get_dataset(step, name='B23',        indices=indices)
        V   = self.get_dataset(step, name='volume',     indices=indices)
        xp  = self.get_dataset(step, name='x_position', indices=indices)
        yp  = self.get_dataset(step, name='y_position', indices=indices)
        zp  = self.get_dataset(step, name='z_position', indices=indices)
        B33 = self._get_B33(B11, B12, B22, B23, V)

        n = len(indices)

        angles = np.empty(n)
        a = np.empty(n)
        b = np.empty(n)
        centres = np.empty((n, 2))

        for i in range(n):
            B = np.array([[B11, B12, B13],
                          [B12, B22, B23],
                          [B13, B23, B33]])

            # calculate inverse of B-matrix: (instead of using np.linalg.inv)
            D, V = np.linalg.eigh(B)
            iD = np.diag(1.0 / D)
            # eigh calculates the eigenvalues in ascending order (c^2, b^2, a^2),
            # hence, the order in iB is (1/a^2, 1/b^2, 1/c^2) assuming a >= b >= c
            iB = np.matmul(np.matmul(V, iD), V.transpose())

            # ellipsoid centre:
            xc = np.array([xp[i], yp[i], zp[i]])

            # create ellipsoid object:
            elid = ellipsoid(xc, iB)

            # calculate ellipes from intersection:
            ell = elid.intersect(pl)

            # get semi-major and semi-minor axes:
            (a[i], b[i]) = ell.semi_axes

            # get ellipse centre:
            (x0, y0) = ell.centre
            centres[i, 0] = x0
            centres[i, 1] = y0

            # get ellipse angle (in degree):
            angles[i] = np.rad2deg(ell.angle)

        return EllipseCollection(widths=2.0 * a,
                                 heights=2.0 * b,
                                 angles=angles,
                                 units='xy',
                                 offsets=centres)

    def _calculated_projected_ellipses(self, step, plane, loc):
        """
        Calculates 2D projections of the ellipsoids onto either
        xy-, xz- or yz-plane.

        Parameters
        ----------
        plane   'xy', 'xz' or 'yz'
        loc     grid point
        """
        if plane not in ['xy', 'xz', 'yz']:
            raise ValueError("Specified plane '", plane "' not available.")

        origin = self.get_origin()
        extent = self.get_box_extent()
        ncells = self.get_box_ncells()

        dx = extent / ncells

        dim = {0: 'x', 1: 'y', 2: 'z'}

        j = 0
        if plane == 'xz':
            j = 1
        elif plane == 'xy':
            j = 2

        # lower and upper position bounds
        lo = origin[j] + dx[j] * (loc - 1)
        hi = origin[j] + dx[j] * (loc + 1)

        # get indices of parcels satisfying lo < pos < hi
        pos = self.get_dataset(step, name=dim[j] + '_position')
        indices = np.where((pos > lo) & (pos < hi))[0]
        pos = None


        B11 = self.get_dataset(step, name='B11',    indices=indices)
        B12 = self.get_dataset(step, name='B12',    indices=indices)
        B13 = self.get_dataset(step, name='B13',    indices=indices)
        B22 = self.get_dataset(step, name='B22',    indices=indices)
        B23 = self.get_dataset(step, name='B23',    indices=indices)
        V   = self.get_dataset(step, name='volume', indices=indices)
        B33 = self._get_B33(B11, B12, B22, B23, V)

        if dim[j] == 'x':
            # yz-plane
            B11 = 1.0 / B11
            B11_proj = B22 - B12 ** 2 * B11
            B12_proj = B23 - B12 * B13 * B11
            B22_proj = B33 - B13 ** 2 * B11
            pos_1 = self.get_dataset(step, name='y_position', indices=indices)
            pos_2 = self.get_dataset(step, name='z_position', indices=indices)
        elif dim[j] == 'y':
            # xz-plane
            B22 = 1.0 / B22
            B11_proj = B11 - B12 ** 2 * B22
            B12_proj = B13 - B12 * B23 * B22
            B22_proj = B33 - B23 ** 2 * B22
            pos_1 = self.get_dataset(step, name='x_position', indices=indices)
            pos_2 = self.get_dataset(step, name='z_position', indices=indices)
        else:
            # xy-plane
            B33 = 1.0 / B33
            B11_proj = B11 - B13 ** 2 * B33
            B12_proj = B12 - B13 * B23 * B33
            B22_proj = B22 - B23 ** 2 * B33
            pos_1 = self.get_dataset(step, name='x_position', indices=indices)
            pos_2 = self.get_dataset(step, name='y_position', indices=indices)

        return pos_1, pos_2, B11_proj, B12_proj, B22_proj



    #def _get_plane(self, data, plane, loc):
        #"""
        #Assumes the data in (x, y, z) ordering.
        #"""
        #if plane not in ['xy', 'xz', 'yz']:
            #raise ValueError("Specified plane '", plane "' not available.")

        #if plane == 'xy':
            #pl_data = data[:, :, loc]
        #elif plane == 'xz':
            #pl_data = data[:, loc, :]
        #else:
            ##plane = 'yz':
            #pl_data = data[loc, :, :]
        #return pl_data
