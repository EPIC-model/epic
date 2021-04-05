AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS = $(top_builddir)/src/parcels/*.o 		\
		$(top_builddir)/src/options.o 		\
		$(top_builddir)/src/writer.o 		\
		$(top_builddir)/src/parameters.o

bin_PROGRAMS = 			\
	ellipse_orientation 	\
	free_slip

ellipse_orientation_SOURCES = ellipse_orientation.f90
free_slip_SOURCES = free_slip.f90
