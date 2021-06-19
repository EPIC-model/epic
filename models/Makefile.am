AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =                                            \
        $(top_builddir)/src/utils/constants.o		\
        $(top_builddir)/src/utils/phys_constants.o	\
        $(top_builddir)/src/utils/options.o		\
        $(top_builddir)/src/utils/parameters.o		\
        $(top_builddir)/src/hdf5/h5_utils.o		\
        $(top_builddir)/src/hdf5/h5_writer.o

modelsdir = $(prefix)/bin

models_PROGRAMS = models

models_SOURCES =		\
	taylorgreen.f90		\
	straka.f90		\
	robert.f90		\
	models.f90

clean-local:
	rm -f *.mod
