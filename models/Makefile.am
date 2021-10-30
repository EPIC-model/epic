AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =                                            \
	$(top_builddir)/src/utils/macros.o		\
        $(top_builddir)/src/utils/constants.o		\
        $(top_builddir)/src/utils/phys_constants.o	\
        $(top_builddir)/src/utils/options.o		\
        $(top_builddir)/src/utils/parameters.o		\
        $(top_builddir)/src/hdf5/h5_utils.o		\
        $(top_builddir)/src/hdf5/h5_writer.o

epic_modelsdir = $(prefix)/bin

epic_models_PROGRAMS = epic-models

epic_models_SOURCES =		\
	taylorgreen.f90		\
	straka.f90		\
	robert.f90		\
	epic-models.f90

clean-local:
	rm -f *.mod
