if ENABLE_UNIT_TESTS

AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS = $(top_builddir)/src/parcels/*.o 		\
		$(top_builddir)/src/options.o 		\
		$(top_builddir)/src/writer.o 		\
		$(top_builddir)/src/parameters.o	\
		$(top_builddir)/src/fields/*.o		\
		$(top_builddir)/src/interpolation/*.o

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 		\
	ellipse_orientation 	\
	free_slip

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

ellipse_orientation_SOURCES = ellipse_orientation.f90
free_slip_SOURCES = free_slip.f90

endif
