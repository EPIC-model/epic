AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =                                            \
        $(top_builddir)/src/utils/constants.o		\
        $(top_builddir)/src/utils/physics.o		\
        $(top_builddir)/src/utils/writer.o

modelsdir = $(prefix)/bin

models_PROGRAMS = models

models_SOURCES =		\
	options.f90		\
	taylorgreen.f90		\
	straka.f90		\
	robert.f90		\
	models.f90

clean-local:
	rm -f *.mod
