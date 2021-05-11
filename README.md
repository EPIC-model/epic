# EPIC - Elliptic Parcel-in-Cell
Elliptic and non-elliptic model for PIC in fluid dynamics

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
* animate-ellipses.py
* plot-ellipses.py
* plot-fields.py
* plot-diagnostics.py
* merge-ellipses.py

Please append the argument `--help` when calling the scripts to get further information. We recommend to install
a separate [Python virtual environment](https://docs.python.org/3/library/venv.html) for EPIC. After activating
the environment all necessary Python packages can be installed via:
```
$ pip install -r requirements.txt
```
The file `requirements.txt` is contained in the root directory of EPIC.

## Performance analysis
When configuring EPIC with `--enable-scalasca`, it is built with the performance tool [Scalasca](https://www.scalasca.org/) and [Score-P](https://www.vi-hps.org/projects/score-p/). Scripts to install Scalasca and Score-P are found in the directory `dependencies`.


