# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/		\
	     -I $(top_builddir)/src/mpi		\
	     -I $(top_builddir)/src/netcdf

lib_LTLIBRARIES = libepic_utils.la
libepic_utils_la_SOURCES = 	\
	datatypes.f90		\
	iomanip.f90		\
	armanip.f90		\
	merge_sort.f90		\
	constants.f90		\
	timer.f90		\
	linalg.f90		\
	jacobi.f90		\
	scherzinger.f90		\
	ape_density.f90		\
	physics.f90

