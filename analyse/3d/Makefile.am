AM_FCFLAGS = 				\
	-I $(top_builddir)/src/3d/	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf	\
	-I $(top_builddir)/src/fft

analysedir = $(prefix)/bin

analyse_PROGRAMS = genspec3d coarsening

genspec3d_SOURCES = genspec3d.f90

genspec3d_LDADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

coarsening_SOURCES = coarsening.f90

coarsening_LDADD =					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod
