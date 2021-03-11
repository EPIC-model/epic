AC_INIT([epic], [1.0], [mf248@st-andrews.ac.uk])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_FC([gfortran])
AC_LANG(Fortran)

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
    Makefile
    src/Makefile
])



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
