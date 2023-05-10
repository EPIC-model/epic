AM_FCFLAGS = -I $(top_builddir)/src/ 	\
	     -I $(top_builddir)/src/mpi

lib_LTLIBRARIES = libepic_netcdf.la
libepic_netcdf_la_SOURCES = 	\
	../config.f90		\
	netcdf_utils.f90	\
	netcdf_writer.f90	\
	netcdf_reader.f90

libepic_netcdf_la_LIBADD = $(top_builddir)/src/mpi/libepic_mpi.la
