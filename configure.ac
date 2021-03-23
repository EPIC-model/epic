AC_INIT([epic], [1.0.0], [mf248@st-andrews.ac.uk], [], [https://github.com/matt-frey/epic])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)


# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

FCFLAGS=""

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
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

AM_PATH_PYTHON([3.0])

#
# numpy
#
AC_MSG_CHECKING([for python numpy library])
`${PYTHON} -c 'import numpy' 2> config.log`
if test "X$?" = "X0" ; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
    AC_MSG_NOTICE([Install numpy with "pip install numpy".])
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
    AC_MSG_NOTICE([Install matplotlib with "pip install matplotlib".])
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
    AC_MSG_NOTICE([Install h5py with "pip install h5py".])
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
    AC_MSG_NOTICE([Install progressbar with "pip install progressbar".])
    exit 1
fi


#######################################################################################
##
## "--enable" flags
##

ENABLE_DEBUG='no'
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug], [compile in debug mode (default=no)])],
              [ENABLE_DEBUG=$enableval])

AM_CONDITIONAL([ENABLE_DEBUG], [test "$ENABLE_DEBUG" = "yes"])

AC_MSG_CHECKING([whether we are compiling in debug mode])
if test "x$ENABLE_DEBUG" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -g -O0 -fbounds-check -fbacktrace"
else
    AC_MSG_RESULT([no])
    FCFLAGS="$FCFLAGS -O3"
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

AC_OUTPUT
