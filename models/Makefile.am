AM_FCFLAGS = 				\
	-I $(top_builddir)/src/2d/	\
	-I $(top_builddir)/src/utils

AM_LDFLAGS =                                            \
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

epic2d_models_LDADD = $(top_builddir)/src/utils/libutils.la

clean-local:
	rm -f *.mod
