if ENABLE_UNIT_TESTS

AM_FCFLAGS = 				\
	-I $(top_builddir)/src/2d/	\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/hdf5

AM_LDFLAGS =						\
	$(top_builddir)/src/2d/parcels/*.o 		\
	$(top_builddir)/src/2d/utils/*.o		\
	$(top_builddir)/src/2d/fields/*.o		\
	$(top_builddir)/src/2d/tri_inversion.o


unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_ellipse_orientation 		\
	test_free_slip				\
	test_fft_1				\
	test_fft_2				\
	test_trilinear				\
	test_laplace_correction			\
	test_gradient_correction		\
	test_parcel_correction		  	\
	test_ellipse_split			\
	test_ellipse_reflection			\
	test_ellipse_merge_1			\
	test_ellipse_merge_2			\
	test_ellipse_merge_3			\
	test_ellipse_merge_4			\
	test_ellipse_merge_5			\
	test_nearest_1				\
	test_nearest_2				\
	test_nearest_3				\
	test_nearest_4				\
	test_tri_inversion			\
	test_parcel_init

test_ellipse_orientation_SOURCES = test_ellipse_orientation.f90
test_ellipse_orientation_LDADD = 			\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_free_slip_SOURCES = test_free_slip.f90
test_free_slip_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_fft_1_SOURCES = test_fft_1.f90
test_fft_1_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_fft_2_SOURCES = test_fft_2.f90
test_fft_2_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_trilinear_SOURCES = test_trilinear.f90
test_trilinear_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_laplace_correction_SOURCES = test_laplace_correction.f90
test_laplace_correction_LDADD = 			\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_gradient_correction_SOURCES = test_gradient_correction.f90
test_gradient_correction_LDADD = 			\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_parcel_correction_SOURCES = test_parcel_correction.f90
test_parcel_correction_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_split_SOURCES = test_ellipse_split.f90
test_ellipse_split_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_reflection_SOURCES = test_ellipse_reflection.f90
test_ellipse_reflection_LDADD = 			\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_merge_1_SOURCES = test_ellipse_merge_1.f90
test_ellipse_merge_1_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_merge_2_SOURCES = test_ellipse_merge_2.f90
test_ellipse_merge_2_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_merge_3_SOURCES = test_ellipse_merge_3.f90
test_ellipse_merge_3_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_merge_4_SOURCES = test_ellipse_merge_4.f90
test_ellipse_merge_4_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_ellipse_merge_5_SOURCES = test_ellipse_merge_5.f90
test_ellipse_merge_5_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_nearest_1_SOURCES = test_nearest_1.f90
test_nearest_1_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_nearest_2_SOURCES = test_nearest_2.f90
test_nearest_2_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_nearest_3_SOURCES = test_nearest_3.f90
test_nearest_3_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_nearest_4_SOURCES = test_nearest_4.f90
test_nearest_4_LDADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_tri_inversion_SOURCES = test_tri_inversion.f90
test_tri_inversion_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

test_parcel_init_SOURCES = test_parcel_init.f90
test_parcel_init_LDADD = 				\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la

if ENABLE_DEBUG
unittests_PROGRAMS += test_ellipse_merge_symmetry
test_ellipse_merge_symmetry_SOURCES = test_ellipse_merge_symmetry.f90
test_ellipse_merge_symmetry_LDADD = 			\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/hdf5/libepic_hdf5.la
endif

endif
