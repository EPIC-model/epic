if ENABLE_UNIT_TESTS
SUBDIRS = 	\
 	.	\
 	2d	\
 	netcdf

if ENABLE_3D
SUBDIRS += 3d
endif

unittestsdir = $(prefix)/unit-tests

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

# 2 November 2021
# https://www.gnu.org/software/libtool/manual/html_node/Using-Automake.html
# https://www.gnu.org/software/libtool/manual/html_node/LT_005fINIT.html
lib_LTLIBRARIES = libunit_test.la
libunit_test_la_SOURCES = unit_test.f90 permute.f90
endif
