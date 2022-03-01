AM_FCFLAGS = 				\
	-I $(top_builddir)/src/2d/	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS = 						\
	$(top_builddir)/src/2d/utils/options.o		\
	$(top_builddir)/src/2d/utils/parameters.o


epic2d_modelsdir = $(prefix)/bin

epic2d_models_PROGRAMS = epic2d-models

epic2d_models_SOURCES =		\
	taylor_green_2d.f90	\
	straka_2d.f90		\
	robert_2d.f90		\
	epic2d-models.f90

epic2d_models_LDADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod
