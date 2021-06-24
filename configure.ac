AC_INIT([epic], [1.0.0], [mf248@st-andrews.ac.uk], [], [https://github.com/matt-frey/epic])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

FCFLAGS="-fdefault-real-8 -fdefault-double-8 -cpp"

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    models/Makefile
    unit-tests/Makefile
    python-scripts/Makefile
    python-scripts/tools/Makefile
    python-scripts/tools/animate/Makefile
])


#######################################################################################
##
## check Python installation
##
## # 23 March 2021
## https://askubuntu.com/questions/29370/how-to-check-if-a-command-succeeded
## https://gitlab.psi.ch/H5hut/src/-/blob/master/configure.ac
##

DISABLE_PYTHON_CHECK='no'
AC_ARG_ENABLE([python-check],
              [AS_HELP_STRING([--disable-python-check], [disable Python checke (default=no)])],
              [DISABLE_PYTHON_CHECK=$disableval])

AM_CONDITIONAL([DISABLE_PYTHON_CHECK], [test "$DISABLE_PYTHON_CHECK" = "yes"])

if test "x$DISABLE_PYTHON_CHECK" = "xno"; then
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
    # h5py
    #
    AC_MSG_CHECKING([for python h5py library])
    `${PYTHON} -c 'import h5py' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install h5py with "conda install -c conda-forge h5py".])
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

#######################################################################################
##
## "--with" flags
##

# 11 March 2021
# https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=<path>], [path to hdf5 library @<:@default=$HDF5_DIR@:>@])],
            [HDF5_DIR=$withval],
            [])

if test -n "${HDF5_DIR}"; then
    FCFLAGS="$FCFLAGS -I$HDF5_DIR/include"
    LDFLAGS="$LDFLAGS -L$HDF5_DIR/lib"
    LIBS="$LIBS -lhdf5_fortran -lhdf5hl_fortran -lhdf5"
fi

# 11 March 2021
# https://github.com/PhysicsofFluids/AFiD/blob/master/configure.ac#L139
AC_MSG_CHECKING([for hdf5 Fortran library])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [
      use hdf5
      implicit none
      integer         :: h5err
      call h5open_f(h5err)])],
  [hdf5_found=yes],
  [hdf5_found=no AC_MSG_ERROR([Cannot find hdf5 Fortran library])])
AC_MSG_RESULT([$hdf5_found])

#######################################################################################
##
## "--enable" flags
##

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
    FCFLAGS="$FCFLAGS -g -O0 -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow"
    FCFLAGS="$FCFLAGS -DENABLE_VERBOSE"
else
    AC_MSG_RESULT([no])
    FCFLAGS="$FCFLAGS -O3 -DNDEBUG"
fi


ENABLE_UNIT_TESTS='no'
AC_ARG_ENABLE([unit-tests],
              [AS_HELP_STRING([--enable-unit-tests], [build with unit tests (default=no)])],
              [ENABLE_UNIT_TESTS=$enableval])

AM_CONDITIONAL([ENABLE_UNIT_TESTS], [test "$ENABLE_UNIT_TESTS" = "yes"])

AC_MSG_CHECKING([whether we build unit tests])
if test "x$ENABLE_UNIT_TESTS" = "xyes"; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
fi

if test "x$ENABLE_UNIT_TESTS" == "xyes" && test "x$ENABLE_DEBUG" == "xno"; then
    AC_MSG_WARN([Volume symmetry unit test disabled since not built in debug mode.])
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
