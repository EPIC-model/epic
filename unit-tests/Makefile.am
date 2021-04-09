if ENABLE_UNIT_TESTS

AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =					\
	$(top_builddir)/src/parcels/*.o 	\
	$(top_builddir)/src/taylorgreen.o	\
	$(top_builddir)/src/utils/*.o		\
	$(top_builddir)/src/fields/*.o		\
	$(top_builddir)/src/interpolation/*.o

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 			\
	test_ellipse_orientation 	\
	test_free_slip			\
	test_fft_1			\
	test_fft_2

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

test_ellipse_orientation_SOURCES = test_ellipse_orientation.f90
test_free_slip_SOURCES = test_free_slip.f90
test_fft_1_SOURCES = test_fft_1.f90
test_fft_2_SOURCES = test_fft_2.f90

endif
