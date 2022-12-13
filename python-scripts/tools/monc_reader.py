from nc_base_reader import nc_base_reader
import numpy as np

class monc_reader(nc_base_reader):

    def get_num_steps(self):
        return self._ncfile.dimensions['time_series_40'].size

    def get_box_extent(self, step=0):
        origin = self.get_box_origin(step)
        x_top = self.get_dataset(step, 'x_top')
        y_top = self.get_dataset(step, 'y_top')
        z_top = self.get_axis(step, 'z')[-1]
        return np.array([x_top, y_top, z_top]) - origin

    def get_box_ncells(self):
        nx = self._ncfile.dimensions['x'].size
        ny = self._ncfile.dimensions['y'].size
        nz = self._ncfile.dimensions['z'].size - 1
        return np.array([nx, ny, nz])

    def get_box_origin(self, step):
        x_bottom = self.get_dataset(step, 'x_bottom')
        y_bottom = self.get_dataset(step, 'y_bottom')
        z_bottom = self.get_axis(step, 'z')[0]
        return np.array([x_bottom, y_bottom, z_bottom])

    def get_axis(self, name, step=0):
        if name == 'x' or name == 'y':
            ng = self._ncfile.dimensions[name].size
            bottom = self.get_dataset(step, name + '_bottom')
            top = self.get_dataset(step, name + '_top')
            return np.linspace(bottom, top, ng+1)
        elif name == 'z':
            return self.get_dataset(step=step, name=name)
        else:
            raise ValueError("No axis called '" + name + "'.")

    def get_dataset(self, step, name, indices=None, copy_periodic=True):
        super().get_dataset(step, name)

        if indices is not None:
            return np.array(self._ncfile.variables[name][step, ...]).squeeze()[indices, ...]
        else:
            fdata = np.array(self._ncfile.variables[name][step, ...]).squeeze()

            if copy_periodic and len(fdata.shape) == 3:
                fdata = self._copy_periodic_layers(fdata)

            return fdata

    def _copy_periodic_layers(self, field):
        nx, ny, nz = field.shape
        field_copy = np.empty((nx+1, ny+1, nz))
        field_copy[0:nx, 0:ny, :] = field.copy()
        field_copy[:, ny, :] = field_copy[:, 0, :]
        field_copy[nx, :, :] = field_copy[0, :, :]
        return field_copy
