AC_INIT([epic], [0.12.4], [mf248@st-andrews.ac.uk], [], [https://github.com/matt-frey/epic])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)
AM_PROG_AR
LT_INIT

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

FCFLAGS="-std=f2018 -fdefault-real-8 -fdefault-double-8 -cpp -fPIC -mcmodel=large -fall-intrinsics"

AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    src/utils/Makefile
    src/mpi/Makefile
    src/netcdf/Makefile
    src/fft/Makefile
    src/2d/Makefile
    src/3d/Makefile
    models/Makefile
    models/2d/Makefile
    models/3d/Makefile
    analyse/Makefile
    analyse/2d/Makefile
    analyse/3d/Makefile
    unit-tests/Makefile
    unit-tests/2d/Makefile
    unit-tests/3d/Makefile
    unit-tests/mpi/Makefile
    unit-tests/netcdf/Makefile
    python-scripts/Makefile
    python-scripts/tools/Makefile
    python-scripts/tools/animate/Makefile
    python-scripts/tools/epic_config.py
    run-script/Makefile
])

#######################################################################################
##
## "--with" flags
##

## define some default search paths
declare default_lib_paths="/lib:/usr/lib:/lib64:/usr/lib64:/usr/local/lib:$LD_LIBRARY_PATH"
declare default_include_paths="/usr/include:/usr/local/lib:$C_INCLUDE_PATH:$CPLUS_INCLUDE_PATH"

# 11 March 2021
# https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<path>], [path to netcdf library @<:@default=$NETCDF_DIR@:>@])],
            [NETCDF_DIR=$withval],
            [])

if test -n "${NETCDF_DIR}"; then
    FCFLAGS="$FCFLAGS -I$NETCDF_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_DIR/lib"
    LIBS="$LIBS -lnetcdff -lnetcdf"
else
    for path in ${default_lib_paths//:/ }; do
        AC_MSG_CHECKING([for netcdf Fortran libraries in $path])
        if test -e "$path/libnetcdff.a"; then
            NETCDF_LIBRARY_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    for path in ${default_include_paths//:/ }; do
        AC_MSG_CHECKING([for netcdf Fortran modules in $path])
        if test -e "$path/netcdf.mod"; then
            NETCDF_INCLUDE_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    if test -n "${NETCDF_LIBRARY_DIR}" && test -n "${NETCDF_INCLUDE_DIR}"; then
        FCFLAGS="$FCFLAGS -I$NETCDF_INCLUDE_DIR"
        LDFLAGS="$LDFLAGS -L$NETCDF_LIBRARY_DIR"
        LIBS="$LIBS -lnetcdff -lnetcdf"
    else
        AC_MSG_ERROR([Cannot find netcdf.
                    Please export the environment variable NETCDF_DIR pointing
                    to the root directory of netcdf.])
    fi
fi

# 11 March 2021
# https://github.com/PhysicsofFluids/AFiD/blob/master/configure.ac#L139
AC_MSG_CHECKING([whether we can compile a netcdf program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    use netcdf
    implicit none])],
[netcdf_found=yes],
[netcdf_found=no AC_MSG_ERROR([Cannot compile a netcdf program])])
AC_MSG_RESULT([$netcdf_found])

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi=<path>], [path to mpi library @<:@default=$MPI_DIR@:>@])],
            [MPI_DIR=$withval],
            [])

if test -n "${MPI_DIR}"; then
    FCFLAGS="$FCFLAGS -I$MPI_DIR/include -I$MPI_DIR/lib"
    LDFLAGS="$LDFLAGS -L$MPI_DIR/lib"
    LIBS="$LIBS -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
else
    for path in ${default_lib_paths//:/ }; do
        AC_MSG_CHECKING([for MPI Fortran libraries in $path])
        if test -e "$path/libmpi_usempif08.a"; then
            MPI_LIBRARY_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    for path in ${default_include_paths//:/ }; do
        AC_MSG_CHECKING([for MPI C header in $path])
        if test -e "$path/mpi.h"; then
            MPI_INCLUDE_DIR="$path"
            AC_MSG_RESULT([yes])
            break
        fi
        AC_MSG_RESULT([no])
    done

    if test -n "${MPI_LIBRARY_DIR}" && test -n "${MPI_INCLUDE_DIR}"; then
        FCFLAGS="$FCFLAGS -I$MPI_INCLUDE_DIR -I$MPI_LIBRARY_DIR"
        LDFLAGS="$LDFLAGS -L$MPI_LIBRARY_DIR"
        LIBS="$LIBS -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
    else
        AC_MSG_ERROR([Cannot find MPI.
                    Please export the environment variable MPI_DIR pointing
                    to the root directory of MPI.])
    fi
fi

AC_MSG_CHECKING([for MPI Fortran library])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [
      use mpi
      implicit none
      integer :: ierr
      call MPI_INIT(ierr)
      call MPI_FINALIZE(ierr)])],
  [mpi_found=yes],
  [mpi_found=no AC_MSG_ERROR([Cannot find mpi Fortran library])])
AC_MSG_RESULT([$mpi_found])

#######################################################################################
##
## "--enable" flags
##

##
## Python
##
## # 23 March 2021
## https://askubuntu.com/questions/29370/how-to-check-if-a-command-succeeded
## https://gitlab.psi.ch/H5hut/src/-/blob/master/configure.ac
##

ENABLE_PYTHON='no'
AC_ARG_ENABLE([python],
              [AS_HELP_STRING([--enable-python], [enable Python (default=no)])],
              [ENABLE_PYTHON=$enableval])

AM_CONDITIONAL([ENABLE_PYTHON], [test "$ENABLE_PYTHON" = "yes"])

if test "x$ENABLE_PYTHON" = "xyes"; then
    AM_PATH_PYTHON([3.0])
    #
    # argparse
    #
    AC_MSG_CHECKING([for python argparse library])
    `${PYTHON} -c 'import argparse' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install argparse with "conda install -c conda-forge argparse".])
        exit 1
    fi

    #
    # numpy
    #
    AC_MSG_CHECKING([for python numpy library])
    `${PYTHON} -c 'import numpy' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install numpy with "conda install -c conda-forge numpy".])
        exit 1
    fi

    #
    # scipy
    #
    AC_MSG_CHECKING([for python scipy library])
    `${PYTHON} -c 'import scipy' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install numpy with "conda install -c conda-forge scipy".])
        exit 1
    fi

    #
    # pandas
    #
    AC_MSG_CHECKING([for python pandas library])
    `${PYTHON} -c 'import pandas' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install numpy with "conda install -c conda-forge pandas".])
        exit 1
    fi

    #
    # matplotlib
    #
    AC_MSG_CHECKING([for python matplotlib library])
    `${PYTHON} -c 'import matplotlib.pyplot as plt' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install matplotlib with "conda install -c conda-forge matplotlib".])
        exit 1
    fi

    #
    # seaborn
    #
    AC_MSG_CHECKING([for python seaborn library])
    `${PYTHON} -c 'import seaborn' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install numpy with "conda install -c conda-forge seaborn".])
        exit 1
    fi

    #
    # netcdf
    #
    AC_MSG_CHECKING([for python netcdf library])
    `${PYTHON} -c 'import netCDF4' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install netcdf with "conda install -c conda-forge netcdf4".])
        exit 1
    fi

    #
    # progressbar
    #
    AC_MSG_CHECKING([for python progressbar library])
    `${PYTHON} -c 'import progressbar' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install progressbar with "conda install -c conda-forge progressbar".])
        exit 1
    fi

    #
    # bokeh
    #
    AC_MSG_CHECKING([for python bokeh library])
    `${PYTHON} -c 'import bokeh' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_NOTICE([Install bokeh with "conda install -c conda-forge bokeh".])
        AC_MSG_NOTICE([Install phantomjs with "conda install -c conda-forge phantomjs".])
        AC_MSG_NOTICE([Install selenium with "conda install -c bokeh selenium".])
        AC_MSG_NOTICE([Install geckodriver with "conda install -c conda-forge firefox geckodriver".])
        exit 1
        AC_MSG_RESULT([no])
    fi
fi

ENABLE_3D='yes'
AC_ARG_ENABLE([3d],
              [AS_HELP_STRING([--enable-3d], [enable 3D model (default=yes)])],
              [ENABLE_3D=$enableval])

AM_CONDITIONAL([ENABLE_3D], [test "$ENABLE_3D" = "yes"])

AC_MSG_CHECKING([whether we are compiling the 3D model])
if test "x$ENABLE_3D" = "xyes"; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
fi

ENABLE_DRY_MODE='no'
AC_ARG_ENABLE([dry-mode],
              [AS_HELP_STRING([--enable-dry-mode], [enable dry mode (default=yes)])],
              [ENABLE_DRY_MODE=$enableval])

AM_CONDITIONAL([ENABLE_DRY_MODE], [test "$ENABLE_DRY_MODE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in dry mode])
if test "x$ENABLE_DRY_MODE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_DRY_MODE"
else
    AC_MSG_RESULT([no])
fi


ENABLE_DIAGNOSE='no'
AC_ARG_ENABLE([diagnose],
              [AS_HELP_STRING([--enable-diagnose], [enable diagnose mode (default=no)])],
              [ENABLE_DIAGNOSE=$enableval])

AM_CONDITIONAL([ENABLE_DIAGNOSE], [test "$ENABLE_DIAGNOSE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in diagnose mode])
if test "x$ENABLE_DIAGNOSE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_DIAGNOSE"
else
    AC_MSG_RESULT([no])
fi

ENABLE_VERBOSE='no'
AC_ARG_ENABLE([verbose],
              [AS_HELP_STRING([--enable-verbose], [enable verbose mode (default=no)])],
              [ENABLE_VERBOSE=$enableval])

AM_CONDITIONAL([ENABLE_VERBOSE], [test "$ENABLE_VERBOSE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in verbose mode])
if test "x$ENABLE_VERBOSE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_VERBOSE"
else
    AC_MSG_RESULT([no])
fi


ENABLE_DEBUG='no'
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug], [compile in debug mode (default=no)])],
              [ENABLE_DEBUG=$enableval])

AM_CONDITIONAL([ENABLE_DEBUG], [test "$ENABLE_DEBUG" = "yes"])

AC_MSG_CHECKING([whether we are compiling in debug mode])
if test "x$ENABLE_DEBUG" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -Wall -Wno-maybe-uninitialized -Werror -g -O0"
    FCFLAGS="$FCFLAGS -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow"
else
    AC_MSG_RESULT([no])
    FCFLAGS="$FCFLAGS -O3 -funroll-all-loops -flto -DNDEBUG"
fi


ENABLE_UNIT_TESTS='no'
AC_ARG_ENABLE([unit-tests],
              [AS_HELP_STRING([--enable-unit-tests], [build with unit tests (default=no)])],
              [ENABLE_UNIT_TESTS=$enableval])

AM_CONDITIONAL([ENABLE_UNIT_TESTS], [test "$ENABLE_UNIT_TESTS" = "yes"])

AC_MSG_CHECKING([whether we build unit tests])
if test "x$ENABLE_UNIT_TESTS" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_UNIT_TESTS"
else
    AC_MSG_RESULT([no])
fi

if test "x$ENABLE_UNIT_TESTS" == "xyes" && test "x$ENABLE_DEBUG" == "xno"; then
    AC_MSG_WARN([Volume symmetry unit test disabled since not built in debug mode.])
fi

ENABLE_OPENMP='no'
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--enable-openmp], [enable OpenMP (default=no)])],
              [ENABLE_OPENMP=$enableval])

AM_CONDITIONAL([ENABLE_OPENMP], [test "$ENABLE_OPENMP" = "yes"])

AC_MSG_CHECKING([whether we are enabling OpenMP])
if test "x$ENABLE_OPENMP" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -fopenmp -DENABLE_OPENMP"
else
    AC_MSG_RESULT([no])
fi

#
# Scalasca and Score-P
#
ENABLE_SCALASCA='no'
AC_ARG_ENABLE([scalasca],
              [AS_HELP_STRING([--enable-scalasca], [compile with scalasca enabled (default=no)])],
              [ENABLE_SCALASCA=$enableval])

AM_CONDITIONAL([ENABLE_SCALASCA], [test "ENABLE_SCALASCA" = "yes"])

AC_MSG_CHECKING([whether we are compiling with scalasca])
if test "x$ENABLE_SCALASCA" = "xyes"; then
    AC_MSG_RESULT([yes])

    # check if scalasca executable is part of $PATH
    AC_CHECK_PROG([scalasca_found],
                   [scalasca],
                   [yes],
                   [no],
                   [$PATH])

    if test "x$scalasca_found" = "xno"; then
        AC_MSG_ERROR([Scalasca not part of PATH environment variable!])
    fi

    # check if scorep (needed by scalasca executable is part of $PATH
    AC_CHECK_PROG([scorep_found],
                   [scorep],
                   [yes],
                   [no],
                   [$PATH])

    if test "x$scorep_found" = "xno"; then
        AC_MSG_ERROR([Score-P not part of PATH environment variable!])
    fi

    FC="scalasca -instrument $FC"
else
    AC_MSG_RESULT([no])
fi

AC_OUTPUT
