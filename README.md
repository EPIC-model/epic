<!--- How to add a license badge found on https://gist.github.com/lukas-h/2a5d00690736b4c3a7ba (1 Feb 2022) --->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5940225.svg)](https://doi.org/10.5281/zenodo.5940225)
[![License](https://img.shields.io/github/license/matt-frey/epic)](https://github.com/matt-frey/epic/blob/main/LICENSE)

# EPIC - Elliptical Parcel-in-Cell
Elliptical PIC model for fluid dynamics

## Dependencies
EPIC has following requirements:
* gfortran
* NetCDF

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

**Note:** You can add the bin directory `$PREFIX/bin` to your `$PATH` and `PYTHONPATH` environment variables with
```
export PATH=$PREFIX/bin:$PATH
export PYTHONPATH=$PREFIX/bin:$PYTHONPATH
```

## Running
In order to run the model, execute
```
$ epic2d --config filename
```
where `filename` is the configuration file. An example of a configuration file is given
[here](examples/taylor_green.config).

## Post-processing
The directory `$PREFIX/bin` contains following Python scripts:
* animate-parcels.py
* plot-parcels.py
* plot-diagnostics.py
* merge-ellipses.py
* subgrid_2d_generalised.py

Please append the argument `--help` when calling the scripts to get further information. We recommend to install
a separate Python virtual environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for EPIC. After installing conda, all necessary Python packages can be installed via
```
$ conda config --add channels conda-forge
$ conda create --name <env> python=3.9.5 --file requirements.txt
```
where `<env>` is the name of the environment. The file `requirements.txt` is contained in the root directory of EPIC.

We further provide Fortran analysis scripts. These are:
* genspec (power spectrum analysis)

## Performance analysis
When configuring EPIC with `--enable-scalasca`, it is built with the performance tool [Scalasca](https://www.scalasca.org/) and [Score-P](https://www.vi-hps.org/projects/score-p/). Scripts to install Scalasca and Score-P are found in the directory `dependencies`.

## How to write an EPIC input field file
EPIC parses a NetCDF file containing all fields to initialise the parcels. You can simply generate the input fields with Python and write them with the provided tools. Below you can find an example where the vorticity field of a Taylor-Green flow is initialised and written to a file.
```Python
#!/usr/bin/env python
#
# Example of writing a field file that can be parsed by EPIC.
#
from tools.nc_fields import nc_fields
import numpy as np

try:
    ncf = nc_fields()

    ncf.open('taylor_green.nc')

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

    vorticity = np.zeros((nz+1, nx))

    # ranges from 0 to nx-1
    for i in range(nx):
        # ranges from 0 to nz
        for j in range(nz+1):
            x = origin[0] + i * dx
            z = origin[1] + j * dz
            vorticity[j, i] = (B * a - A * b) * np.cos(a * x + d) * np.cos(b * z + e)

    # write all provided fields
    ncf.add_field('vorticity', vorticity, unit='1/s')

    ncf.add_box(origin, extent, [nx, nz])

    ncf.close()

except Exception as err:
    print(err)
```
