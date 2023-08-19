AM_FCFLAGS = 				\
	-I $(top_builddir)/src/2d/	\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =						\
	$(top_builddir)/src/2d/parcels/*.o 		\
	$(top_builddir)/src/2d/utils/*.o		\
	$(top_builddir)/src/2d/fields/*.o		\
	$(top_builddir)/src/2d/tri_inversion.o

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_ellipse_orientation 		\
	test_free_slip				\
	test_fft_1				\
	test_fft_2				\
	test_bilinear				\
	test_laplace_correction			\
	test_gradient_correction		\
	test_parcel_correction		  	\
	test_ellipse_split			\
	test_ellipse_reflection			\
	test_ellipse_merge_1			\
	test_ellipse_merge_2			\
	test_ellipse_merge_3			\
	test_ellipse_merge_4			\
	test_nearest_1				\
	test_nearest_2				\
	test_nearest_3				\
	test_nearest_4				\
	test_tri_inversion			\
	test_parcel_init

test_ellipse_orientation_SOURCES = test_ellipse_orientation.f90
test_ellipse_orientation_LDADD = libcombi.la

test_free_slip_SOURCES = test_free_slip.f90
test_free_slip_LDADD = libcombi.la

test_fft_1_SOURCES = test_fft_1.f90
test_fft_1_LDADD = libcombi.la

test_fft_2_SOURCES = test_fft_2.f90
test_fft_2_LDADD = libcombi.la

test_bilinear_SOURCES = test_bilinear.f90
test_bilinear_LDADD = libcombi.la

test_laplace_correction_SOURCES = test_laplace_correction.f90
test_laplace_correction_LDADD = libcombi.la

test_gradient_correction_SOURCES = test_gradient_correction.f90
test_gradient_correction_LDADD = libcombi.la

test_parcel_correction_SOURCES = test_parcel_correction.f90
test_parcel_correction_LDADD = libcombi.la

test_ellipse_split_SOURCES = test_ellipse_split.f90
test_ellipse_split_LDADD = libcombi.la

test_ellipse_reflection_SOURCES = test_ellipse_reflection.f90
test_ellipse_reflection_LDADD = libcombi.la

test_ellipse_merge_1_SOURCES = test_ellipse_merge_1.f90
test_ellipse_merge_1_LDADD = libcombi.la

test_ellipse_merge_2_SOURCES = test_ellipse_merge_2.f90
test_ellipse_merge_2_LDADD = libcombi.la

test_ellipse_merge_3_SOURCES = test_ellipse_merge_3.f90
test_ellipse_merge_3_LDADD = libcombi.la

test_ellipse_merge_4_SOURCES = test_ellipse_merge_4.f90
test_ellipse_merge_4_LDADD = libcombi.la

test_nearest_1_SOURCES = test_nearest_1.f90
test_nearest_1_LDADD = libcombi.la
test_nearest_2_SOURCES = test_nearest_2.f90
test_nearest_2_LDADD = libcombi.la

test_nearest_3_SOURCES = test_nearest_3.f90
test_nearest_3_LDADD = libcombi.la

test_nearest_4_SOURCES = test_nearest_4.f90
test_nearest_4_LDADD = libcombi.la

test_tri_inversion_SOURCES = test_tri_inversion.f90
test_tri_inversion_LDADD = libcombi.la
test_parcel_init_SOURCES = test_parcel_init.f90
test_parcel_init_LDADD = libcombi.la

