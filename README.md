# EPIC - Elliptic Parcel-in-Cell
Elliptic and non-elliptic model for PIC in fluid dynamics

## Dependencies
EPIC has following requirements:
* gfortran
* hdf5

A script to install hdf5 is found in the subdirectory `dependencies`. If you do not install hdf5 to the system
location, replace `../configure` below by
```
$ ../configure --with-hdf5=$PREFIX
```
where `$PREFIX` is the root directory of your hdf5 installation.

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
* merge-ellipses.py

Please append the argument `--help` when calling the scripts to get further information.


## Performance analysis
When configuring EPIC with `--enable-scalasca`, it is built with the performance tool [Scalasca](https://www.scalasca.org/) and [Score-P](https://www.vi-hps.org/projects/score-p/). Scripts to install Scalasca and Score-P are found in the directory `dependencies`.


