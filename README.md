# EPIC - Elliptical Parcel-in-Cell
Elliptical PIC model for fluid dynamics

## Dependencies
EPIC has following requirements:
* gfortran
* zlib
* hdf5

The scripts to install hdf5 and zlib are found in the subdirectory `dependencies`. If you do not install hdf5 to
the system location, replace `../configure --prefix=$HDF5` below by
```
$ ../configure --prefix=$PREFIX --with-hdf5=$HDF5
```
where `$HDF5` is the root directory of your hdf5 installation. To install hdf5 with zlib, you need to export
the variable `ZLIB_DIR` to the root install directory.

## Compile
In the following `$PREFIX` denotes the installation directory of EPIC.
Type
```
$ autoreconf --install
$ mkdir build
$ cd build
$ ../configure --prefix=$PREFIX
$ make
$ make install
```

**Note:** You can add the bin directory `$PREFIX/bin` to your `$PATH` environment variable with
```
export PATH=$PREFIX/bin:$PATH
```

## Running
In order to run the model, execute
```
$ epic --config filename
```
where `filename` is the configuration file.

## Post-processing
The directory `$PREFIX/bin` contains following scripts:
* animate-parcels.py
* plot-parcels.py
* plot-diagnostics.py
* merge-ellipses.py

Please append the argument `--help` when calling the scripts to get further information. We recommend to install
a separate Python virtual environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for EPIC. After installing conda, all necessary Python packages can be installed via
```
$ conda create --name <env> --file requirements.txt
```
where `<env>` is the name of the environment. The file `requirements.txt` is contained in the root directory of EPIC.

## Performance analysis
When configuring EPIC with `--enable-scalasca`, it is built with the performance tool [Scalasca](https://www.scalasca.org/) and [Score-P](https://www.vi-hps.org/projects/score-p/). Scripts to install Scalasca and Score-P are found in the directory `dependencies`.

## How to write an EPIC input field file
EPIC parses a HDF5 file containing all fields to initialise the parcels. You can simply generate the input fields with Python and write them with the provided tools. Below you can find an example where the vorticity field of a Taylor-Green flow is initialised and written to a file.
```Python
#!/usr/bin/env python
#
# Example of writing a field file that can be parsed by EPIC.
#
from h5_fields import h5fields
import numpy as np

try:
    h5f = h5fields()

    h5f.open('taylor_green.hdf5')

    # velocity field:
    # u(x, z) = A * cos(ax + d) * sin(bz + e)
    # w(x, z) = B * sin(ax + d) * cos(bz + e)

    # vorticity:
    # zeta = (B * a - A * b) * cos(ax + d) * cos(bz + e)

    # amplitudes
    A = 0.5
    B = -1.0

    # frequencies
    a = 2.0
    b = 1.0

    # phases
    d = 0.5 * np.pi
    e = 0.0

    # number of cells
    nx = 32
    nz = 32

    # domain origin
    origin = (-0.5 * np.pi, -0.5 * np.pi)

    # domain extent
    extent = (np.pi, np.pi)

    # mesh spacings
    dx = extent[0] / nx
    dz = extent[1] / nz

    vorticity = np.zeros((nx, nz+1))

    # ranges from 0 to nx-1
    for i in range(nx):
        # ranges from 0 to nz
        for j in range(nz+1):
            x = origin[0] + i * dx
            z = origin[0] + j * dz
            vorticity[i, j] = (B * a - A * b) * np.cos(a * x + d) * np.cos(b * z + e)

    fields = {'vorticity': vorticity}

    # write all provided fields
    h5f.add_fields(origin, extent, fields)

    h5f.close()

except Exception as err:
    print(err)
```


