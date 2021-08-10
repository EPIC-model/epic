if ENABLE_UNIT_TESTS

AM_FCFLAGS = -I $(top_builddir)/src/ -Wno-uninitialized

AM_LDFLAGS =						\
	$(top_builddir)/src/parcels/*.o 		\
	$(top_builddir)/src/utils/*.o			\
	$(top_builddir)/src/hdf5/*.o			\
	$(top_builddir)/src/fields/*.o			\
	$(top_builddir)/src/tri_inversion.o

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
	test_tri_inversion			\
	test_parcel_init

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

test_ellipse_orientation_SOURCES = unit_test.f90 test_ellipse_orientation.f90
test_free_slip_SOURCES = test_free_slip.f90 unit_test.f90
test_fft_1_SOURCES = test_fft_1.f90 unit_test.f90
test_fft_2_SOURCES = test_fft_2.f90 unit_test.f90
test_trilinear_SOURCES = test_trilinear.f90 unit_test.f90
test_laplace_correction_SOURCES = test_laplace_correction.f90 unit_test.f90
test_gradient_correction_SOURCES = test_gradient_correction.f90 unit_test.f90
test_parcel_correction_SOURCES = test_parcel_correction.f90 unit_test.f90
test_ellipse_split_SOURCES = test_ellipse_split.f90 unit_test.f90
test_ellipse_reflection_SOURCES = test_ellipse_reflection.f90 unit_test.f90
test_ellipse_merge_1_SOURCES = test_ellipse_merge_1.f90 unit_test.f90
test_ellipse_merge_2_SOURCES = test_ellipse_merge_2.f90 unit_test.f90
test_ellipse_merge_3_SOURCES = test_ellipse_merge_3.f90 unit_test.f90
test_ellipse_merge_4_SOURCES = test_ellipse_merge_4.f90 unit_test.f90
test_ellipse_merge_5_SOURCES = test_ellipse_merge_5.f90 unit_test.f90
test_nearest_1_SOURCES = test_nearest_1.f90 unit_test.f90
test_nearest_2_SOURCES = test_nearest_2.f90 unit_test.f90
test_tri_inversion_SOURCES = test_tri_inversion.f90 unit_test.f90
test_parcel_init_SOURCES = test_parcel_init.f90 unit_test.f90



if ENABLE_DEBUG
unittests_PROGRAMS += test_ellipse_merge_symmetry
test_ellipse_merge_symmetry_SOURCES = test_ellipse_merge_symmetry.f90 unit_test.f90
endif

endif
