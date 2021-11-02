AM_FCFLAGS = -I $(top_builddir)/src/2d/

AM_LDFLAGS =                                            \
        $(top_builddir)/src/2d/utils/constants.o	\
        $(top_builddir)/src/2d/utils/phys_constants.o	\
        $(top_builddir)/src/2d/utils/options.o		\
        $(top_builddir)/src/2d/utils/parameters.o	\
        $(top_builddir)/src/2d/hdf5/h5_utils.o		\
        $(top_builddir)/src/2d/hdf5/h5_writer.o

epic2d_modelsdir = $(prefix)/bin

epic2d_models_PROGRAMS = epic2d-models

epic2d_models_SOURCES =		\
	taylorgreen.f90		\
	straka.f90		\
	robert.f90		\
	epic2d-models.f90

clean-local:
	rm -f *.mod
