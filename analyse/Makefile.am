AM_FCFLAGS = -I $(top_builddir)/src/2d/

AM_LDFLAGS =                                            \
        $(top_builddir)/src/2d/utils/constants.o	\
        $(top_builddir)/src/2d/utils/stafft.o		\
        $(top_builddir)/src/2d/utils/sta2dfft.o		\
        $(top_builddir)/src/2d/utils/options.o		\
        $(top_builddir)/src/2d/utils/parameters.o      	\
        $(top_builddir)/src/2d/hdf5/h5_utils.o		\
        $(top_builddir)/src/2d/hdf5/h5_reader.o		\
        $(top_builddir)/src/2d/hdf5/h5_writer.o

analysedir = $(prefix)/bin

analyse_PROGRAMS = genspec

genspec_SOURCES = genspec.f90

clean-local:
	rm -f *.mod
