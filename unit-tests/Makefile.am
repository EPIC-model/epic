SUBDIRS = 	\
	.	\
	2d	\
	3d

AM_FCFLAGS =				\
	-I $(top_builddir)/src/2d 	\
	-I $(top_builddir)/src/utils

# AM_LDFLAGS =                                            \
        $(top_builddir)/src/2d/utils/*.o



unittestsdir = $(prefix)/unit-tests

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

# 2 November 2021
# https://www.gnu.org/software/libtool/manual/html_node/Using-Automake.html
# https://www.gnu.org/software/libtool/manual/html_node/LT_005fINIT.html
lib_LTLIBRARIES = libunit_test.la
libunit_test_la_SOURCES = unit_test.f90 permute.f90
