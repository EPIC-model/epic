AC_INIT([epic], [0.14.3], [mf248@st-andrews.ac.uk], [], [https://github.com/EPIC-model/epic])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

: ${CFLAGS=""}
: ${FCFLAGS=""}
: ${CXXFLAGS=""}

AC_PROG_FC([gfortran])
AC_LANG(Fortran)
AM_PROG_AR
LT_INIT

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)
AC_FC_PP_SRCEXT(f90)

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
    mpi-tests/Makefile
    python-scripts/Makefile
    python-scripts/tools/Makefile
    python-scripts/tools/animate/Makefile
    python-scripts/tools/epic_config.py
    run-script/Makefile
])

# 30 May 2023
# https://www.gnu.org/software/autoconf-archive/ax_compiler_vendor.html
# https://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
AX_COMPILER_VENDOR

if test "$ax_cv_fc_compiler_vendor" = "unknown"; then
    AC_MSG_ERROR([Could not deduce compiler vendor.])
fi

AX_CHECK_COMPILE_FLAG([-cpp], [FCFLAGS="$FCFLAGS -cpp"])
AX_CHECK_COMPILE_FLAG([-mcmodel=large], [FCFLAGS="$FCFLAGS -mcmodel=large"])

# gfortran compiler flags
if test "$ax_cv_fc_compiler_vendor" = "gnu"; then
   AX_CHECK_COMPILE_FLAG([-std=f2018], [FCFLAGS="$FCFLAGS -std=f2018"])
   AX_CHECK_COMPILE_FLAG([-fall-intrinsics], [FCFLAGS="$FCFLAGS -fall-intrinsics"])
   AX_CHECK_COMPILE_FLAG([-fdefault-real-8], [FCFLAGS="$FCFLAGS -fdefault-real-8"])
   AX_CHECK_COMPILE_FLAG([-fdefault-double-8], [FCFLAGS="$FCFLAGS -fdefault-double-8"])
   AX_CHECK_COMPILE_FLAG([-Wno-unused-dummy-argument], [FCFLAGS="$FCFLAGS -Wno-unused-dummy-argument"])
fi

# intel compiler flags
AX_CHECK_COMPILE_FLAG([-std18], [FCFLAGS="$FCFLAGS -std18"])
AX_CHECK_COMPILE_FLAG([-fp-model=source], [FCFLAGS="$FCFLAGS -fp-model=source"])
AX_CHECK_COMPILE_FLAG([-real-size 8], [FCFLAGS="$FCFLAGS -real-size 8"])
AX_CHECK_COMPILE_FLAG([-double-size 64], [FCFLAGS="$FCFLAGS -double-size 64"])


#######################################################################################
##
## "--with" flags
##


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF C

## 11 March 2021
## https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([netcdf-c],
            [AS_HELP_STRING([--with-netcdf-c=<path>], [path to NetCDF C library @<:@default=$NETCDF_C_DIR@:>@])],
            [NETCDF_C_DIR=$withval],
            [])

if test -n "${NETCDF_C_DIR}"; then
    CFLAGS="$CFLAGS -I$NETCDF_C_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_C_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_C_DIR environment variable.])
fi

AC_LANG_PUSH([C])
AC_SEARCH_LIBS([nc_create], [netcdf], [], [], [])
AC_MSG_CHECKING([whether we can compile a NetCDF C program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([#include <netcdf.h>], [
    int ncid;
    nc_create("dummy.nc", NC_CLOBBER, &ncid);
    ])],
[netcdf_c_found=yes],
[netcdf_c_found=no])
AC_MSG_RESULT([$netcdf_c_found])
AC_LANG_POP([C])

if test "x$netcdf_c_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF C program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF Fortran

AC_ARG_WITH([netcdf-fortran],
            [AS_HELP_STRING([--with-netcdf-fortran=<path>],
            [path to NetCDF Fortran library @<:@default=$NETCDF_FORTRAN_DIR@:>@])],
            [NETCDF_FORTRAN_DIR=$withval],
            [])

if test -n "${NETCDF_FORTRAN_DIR}"; then
    FCFLAGS="$FCFLAGS -I$NETCDF_FORTRAN_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_FORTRAN_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_FORTRAN_DIR environment variable.])
fi


# 15 April 2023
#http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_prog_fc_mpi.m4
# We cannot use AC_SEARCH_LIBS
for lib in NONE netcdff; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for function nf90_open])
    else
    AC_MSG_CHECKING([for function nf90_open in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
           use netcdf
            implicit none
            integer :: ncerr, ncid
            ncerr = nf90_open("dummy.nc", NF90_NOWRITE, ncid)])],
        [netcdf_found=yes],
        [netcdf_found=no])
    AC_MSG_RESULT($netcdf_found)
    if test "x$netcdf_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done

AC_MSG_CHECKING([whether we can compile a NetCDF Fortran program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    use netcdf
    implicit none])],
[netcdf_fortran_found=yes],
[netcdf_fortran_found=no])
AC_MSG_RESULT([$netcdf_fortran_found])

if test "x$netcdf_fortran_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF Fortran program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MPI

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi=<path>], [path to mpi library @<:@default=$MPI_DIR@:>@])],
            [MPI_DIR=$withval],
            [])

if test -n "${MPI_DIR}"; then
    # mod file may also be in $MPI_DIR/lib
    FCFLAGS="$FCFLAGS -I$MPI_DIR/include -I$MPI_DIR/lib"
    LDFLAGS="$LDFLAGS -L$MPI_DIR/lib"
else
    AC_MSG_ERROR([No MPI_DIR environment variable.])
fi

AC_SEARCH_LIBS([MPI_Init],
               [mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90],
               [],
               [])

AC_SEARCH_LIBS([MPI_Win_lock],
               [mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90],
               [],
               [])

for lib in NONE mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for module mpi_f08])
    else
    AC_MSG_CHECKING([for module mpi_f08 in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
        [mpi_found=yes],
        [mpi_found=no])
    AC_MSG_RESULT([$mpi_found])

    if test "x$mpi_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done


AC_MSG_CHECKING([for MPI Fortran library])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
    [mpi_found=yes],
    [mpi_found=no])
AC_MSG_RESULT([$mpi_found])

if test "x$mpi_found" = "xno"; then
    AC_MSG_ERROR([Cannot find MPI Fortran library])
fi

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

ENABLE_PYTHON_CHECK='no'
AC_ARG_ENABLE([python],
              [AS_HELP_STRING([--enable-python-check], [enable Python check (default=no)])],
              [ENABLE_PYTHON_CHECK=$enableval])

AM_CONDITIONAL([ENABLE_PYTHON_CHECK], [test "$ENABLE_PYTHON_CHECK" = "yes"])

if test "x$ENABLE_PYTHON_CHECK" = "xyes"; then
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

ENABLE_BUOYANCY_PERTURBATION_MODE='no'
AC_ARG_ENABLE([buoyancy-perturbation-mode],
              [AS_HELP_STRING([--enable-buoyancy-perturbation-mode],
              [enable buoyancy perturbation mode (default=no)])],
              [ENABLE_BUOYANCY_PERTURBATION_MODE=$enableval])

AM_CONDITIONAL([ENABLE_BUOYANCY_PERTURBATION_MODE], [test "$ENABLE_BUOYANCY_PERTURBATION_MODE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in buoyancy perturbation mode])
if test "x$ENABLE_BUOYANCY_PERTURBATION_MODE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_BUOYANCY_PERTURBATION_MODE"
else
    AC_MSG_RESULT([no])
fi

ENABLE_DRY_MODE='no'
AC_ARG_ENABLE([dry-mode],
              [AS_HELP_STRING([--enable-dry-mode], [enable dry mode (default=no)])],
              [ENABLE_DRY_MODE=$enableval])

AM_CONDITIONAL([ENABLE_DRY_MODE], [test "$ENABLE_DRY_MODE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in dry mode])
if test "x$ENABLE_DRY_MODE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_DRY_MODE"
else
    AC_MSG_RESULT([no])
fi

ENABLE_LABELS='no'
AC_ARG_ENABLE([labels],
              [AS_HELP_STRING([--enable-labels], [enable labels (default=no)])],
              [ENABLE_LABELS=$enableval])

AM_CONDITIONAL([ENABLE_LABELS], [test "$ENABLE_LABELS" = "yes"])

AC_MSG_CHECKING([whether we are compiling with labels])
if test "x$ENABLE_LABELS" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_LABELS"
else
    AC_MSG_RESULT([no])
fi

ENABLE_P2G_1POINT='no'
AC_ARG_ENABLE([p2g-1point],
              [AS_HELP_STRING([--enable-p2g-1point], [enable one-point par2grid (default=no)])],
              [ENABLE_P2G_1POINT=$enableval])

AM_CONDITIONAL([ENABLE_P2G_1POINT], [test "$ENABLE_P2G_1POINT" = "yes"])

AC_MSG_CHECKING([whether we are compiling with 1 point in par2grid])
if test "x$ENABLE_P2G_1POINT" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_P2G_1POINT"
else
    AC_MSG_RESULT([no])
fi

ENABLE_G2P_1POINT='no'
AC_ARG_ENABLE([g2p-1point],
              [AS_HELP_STRING([--enable-g2p-1point], [enable one-point grid2par (default=no)])],
              [ENABLE_G2P_1POINT=$enableval])

AM_CONDITIONAL([ENABLE_G2P_1POINT], [test "$ENABLE_G2P_1POINT" = "yes"])

AC_MSG_CHECKING([whether we are compiling with 1 point in grid2par])
if test "x$ENABLE_G2P_1POINT" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_G2P_1POINT"
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
    AX_CHECK_COMPILE_FLAG([-g], [FCFLAGS="$FCFLAGS -g"])
    AX_CHECK_COMPILE_FLAG([-O0], [FCFLAGS="$FCFLAGS -O0"])

    # gfortran compiler flags
    AX_CHECK_COMPILE_FLAG([-Wall], [FCFLAGS="$FCFLAGS -Wall"])
    AX_CHECK_COMPILE_FLAG([-Wuninitialized], [FCFLAGS="$FCFLAGS -Wuninitialized"])
    AX_CHECK_COMPILE_FLAG([-Wmaybe-uninitialized], [FCFLAGS="$FCFLAGS -Wmaybe-uninitialized"])
    AX_CHECK_COMPILE_FLAG([-Werror], [FCFLAGS="$FCFLAGS -Werror"])
    AX_CHECK_COMPILE_FLAG([-fcheck=all], [FCFLAGS="$FCFLAGS -fcheck=all"])
    AX_CHECK_COMPILE_FLAG([-fbounds-check], [FCFLAGS="$FCFLAGS -fbounds-check"])
    AX_CHECK_COMPILE_FLAG([-fbacktrace], [FCFLAGS="$FCFLAGS -fbacktrace"])
    AX_CHECK_COMPILE_FLAG([-ffpe-trap=denormal,invalid,zero,overflow,underflow],
                          [FCFLAGS="$FCFLAGS -ffpe-trap=denormal,invalid,zero,overflow,underflow"])

    # intel compiler flags
    AX_CHECK_COMPILE_FLAG([-warn all], [FCFLAGS="$FCFLAGS -Wall"])
    AX_CHECK_COMPILE_FLAG([-warn error], [FCFLAGS="$FCFLAGS -warn error"])
    AX_CHECK_COMPILE_FLAG([-debug full], [FCFLAGS="$FCFLAGS -debug full"])
else
    AC_MSG_RESULT([no])
    AX_CHECK_COMPILE_FLAG([-O3], [FCFLAGS="$FCFLAGS -O3"])
    AX_CHECK_COMPILE_FLAG([-DNDEBUG], [FCFLAGS="$FCFLAGS -DNDEBUG"])

    # gfortran compiler flags
    AX_CHECK_COMPILE_FLAG([-funroll-all-loops], [FCFLAGS="$FCFLAGS -funroll-all-loops"])
    AX_CHECK_COMPILE_FLAG([-flto], [FCFLAGS="$FCFLAGS -flto"])

    # intel compiler flags
    AX_CHECK_COMPILE_FLAG([-funroll-loops], [FCFLAGS="$FCFLAGS -funroll-loops"])
    AX_CHECK_COMPILE_FLAG([-ipo -c], [FCFLAGS="$FCFLAGS -ipo -c"])
    AX_CHECK_COMPILE_FLAG([-ffat-lto-objects], [FCFLAGS="$FCFLAGS -ffat-lto-objects"])
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
    # gfortran compiler flags
    AX_CHECK_COMPILE_FLAG([-fopenmp], [FCFLAGS="$FCFLAGS -fopenmp -DENABLE_OPENMP"])

    # intel compiler flags
    AX_CHECK_COMPILE_FLAG([-qopenmp], [FCFLAGS="$FCFLAGS -qopenmp -DENABLE_OPENMP"])
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
