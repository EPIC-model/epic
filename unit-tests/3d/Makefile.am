if ENABLE_UNIT_TESTS

AM_FCFLAGS =				\
	-I $(top_builddir)/src/3d/ 	\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/hdf5

AM_LDFLAGS =						\
	$(top_builddir)/src/3d/parcels/*.o 		\
	$(top_builddir)/src/3d/utils/*.o		\
	$(top_builddir)/src/3d/fields/*.o		\
	$(top_builddir)/src/3d/tri_inversion.o

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_jacobi_1				\
	test_jacobi_2				\
	test_jacobi_3				\
	test_ellipsoid_split \
	test_parcel_init_3d

dataroot_DATA =		\
	numpy_eigh.py 	\
	V.asc	   	\
	D.asc	   	\
	B.asc

test_jacobi_1_SOURCES = test_jacobi_1.f90
test_jacobi_1_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_jacobi_2_SOURCES = test_jacobi_2.f90
test_jacobi_2_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_jacobi_3_SOURCES = test_jacobi_3.f90
test_jacobi_3_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipsoid_split_SOURCES = test_ellipsoid_split.f90
test_ellipsoid_split_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_parcel_init_3d_SOURCES = test_parcel_init_3d.f90
test_parcel_init_3d_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

endif
