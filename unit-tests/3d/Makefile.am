if ENABLE_UNIT_TESTS

AM_FCFLAGS = -I $(top_builddir)/src/ -I $(top_builddir)/unit-tests

AM_LDFLAGS =						\
	$(top_builddir)/src/parcels/*.o 		\
	$(top_builddir)/src/utils/*.o			\
	$(top_builddir)/src/hdf5/*.o			\
	$(top_builddir)/src/fields/*.o			\
	$(top_builddir)/src/tri_inversion.o

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_jacobi_1				\
	test_jacobi_2				\
	test_jacobi_3

dataroot_DATA =		\
	numpy_eigh.py 	\
	V.asc	   	\
	D.asc	   	\
	B.asc

test_jacobi_1_SOURCES = test_jacobi_1.f90
test_jacobi_1_LDADD = $(top_builddir)/unit-tests/libunit_test.la

test_jacobi_2_SOURCES = test_jacobi_2.f90
test_jacobi_2_LDADD = $(top_builddir)/unit-tests/libunit_test.la

test_jacobi_3_SOURCES = test_jacobi_3.f90
test_jacobi_3_LDADD = $(top_builddir)/unit-tests/libunit_test.la


endif
