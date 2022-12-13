from nc_base_reader import nc_base_reader
import numpy as np

class pmpic_reader(nc_base_reader):

    @property
    def _has_time_dimension(self):
        return ('time' in self._ncfile.dimensions)

    def get_num_steps(self):
        if self._has_time_dimension:
            return self._ncfile.dimensions['time'].size
        return 1

    def get_box_extent(self):
        origin = self.get_box_origin()
        x_top = self.get_axis('x')[-1]
        y_top = self.get_axis('y')[-1]
        z_top = self.get_axis('z')[-1]
        return np.array([x_top, y_top, z_top]) - origin

    def get_box_ncells(self):
        nx = self._ncfile.dimensions['xp'].size
        ny = self._ncfile.dimensions['yp'].size
        nz = self._ncfile.dimensions['zp'].size - 1
        return np.array([nx, ny, nz])

    def get_box_origin(self):
        x_origin = self.get_axis('x')[0]
        y_origin = self.get_axis('y')[0]
        z_origin = self.get_axis('z')[0]
        return np.array([x_origin, y_origin, z_origin])

    def get_axis(self, name):
        if not name in ['x', 'y', 'z']:
            raise ValueError("No axis called '" + name + "'.")
        n = name + 'p'
        axis = self.get_all(n)
        if name == 'x' or name == 'y':
            delta = axis[1] - axis[0]
            new_axis = np.empty(axis.size + 1)
            new_axis[0:-1] = axis[:]
            new_axis[-1] = axis[-1] + delta
            axis = new_axis
        return axis

    def get_dataset(self, step, name):
        super().get_dataset(step, name)

        if self._has_time_dimension:
            fdata = np.array(self._ncfile.variables[name][step, ...])
        else:
            fdata = np.array(self._ncfile.variables[name])
        fdata = self._copy_periodic_layers(fdata)
        return fdata

    def _copy_periodic_layers(self, field):
        nx, ny, nz = field.shape
        field_copy = np.empty((nx+1, ny+1, nz))
        field_copy[0:nx, 0:ny, :] = field.copy()
        field_copy[:, ny, :] = field_copy[:, 0, :]
        field_copy[nx, :, :] = field_copy[0, :, :]
        return field_copy
