# EPIC - Elliptic Parcel-in-Cell
Elliptic and non-elliptic model for PIC in fluid dynamics


## Building HDF5 from scratch
Tested with
zlib-1.2.11 (https://zlib.net/)
hdf5-1.10.5 (https://support.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.10/hdf5-1.10.5/src/)


Example compilation command showing some possible dependencies (use sudo if needed):
For zlib
```
$ ./configure --prefix=$HOME/dependencies/zlib CC=mpicc
$ sudo make install
```

For HDF5 with fortran enabled
```
$ ./configure --with-zlib=$HOME/dependencies/zlib --enable-parallel --enable-fortran --prefix=$HOME/dependencies/hdf5 CC='/usr/lib64/openmpi/bin/mpicc' CXX='/usr/lib64/openmpi/bin/mpicxx' FC='/usr/lib64/openmpi/bin/mpif90'
$ sudo make install
```

In this case, also replace ../configure below by
```
$ ../configure --with-hdf5="$HOME/dependencies/hdf5"
```

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
$ ./epic --config filename
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


