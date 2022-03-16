AM_FCFLAGS = 				\
	-I $(top_builddir)/src/2d/	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =                                            \
        $(top_builddir)/src/2d/utils/stafft.o		\
        $(top_builddir)/src/2d/utils/sta2dfft.o

analysedir = $(prefix)/bin

analyse_PROGRAMS = genspec

genspec_SOURCES = genspec.f90

genspec_LDADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod
