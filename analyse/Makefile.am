AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =                                            \
	$(top_builddir)/src/utils/macros.o		\
        $(top_builddir)/src/utils/constants.o		\
        $(top_builddir)/src/utils/stafft.o		\
        $(top_builddir)/src/utils/sta2dfft.o		\
        $(top_builddir)/src/utils/options.o		\
        $(top_builddir)/src/utils/parameters.o          \
        $(top_builddir)/src/hdf5/h5_utils.o		\
        $(top_builddir)/src/hdf5/h5_reader.o		\
        $(top_builddir)/src/hdf5/h5_writer.o

analysedir = $(prefix)/bin

analyse_PROGRAMS = genspec

genspec_SOURCES = genspec.f90

clean-local:
	rm -f *.mod
