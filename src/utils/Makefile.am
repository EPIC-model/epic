# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/utils	\
	     -I $(top_builddir)/src/netcdf

lib_LTLIBRARIES = libepic_utils.la
libepic_utils_la_SOURCES = 	\
	merge_sort.f90		\
	constants.f90		\
	physical_constants.f90	\
	timer.f90		\
	jacobi.f90		\
	iomanip.f90		\
	config.f90		\
	physical_parameters.f90

