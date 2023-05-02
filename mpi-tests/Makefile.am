if ENABLE_3D
mpitestsdir = $(prefix)/mpi-tests

AM_FCFLAGS =				\
	-I $(top_builddir)/src/3d/ 	\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =						\
	$(top_builddir)/src/3d/parcels/*.o 		\
	$(top_builddir)/src/3d/utils/*.o		\
	$(top_builddir)/src/3d/fields/*.o		\
	$(top_builddir)/src/3d/inversion/*.o		\
	$(top_builddir)/src/3d/stepper/*.o

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/mpi/libepic_mpi.la		\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/fft/libepic_mpi_fft.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

mpitests_PROGRAMS =			\
	test_parcel_moving_equidistant	\
	test_parcel_moving_random	\
	test_parcel_split_random	\
	test_parcel_merge_random	\
	test_parcel_split_merge		\
	test_merging_parcels		\
	test_merging_parcels_serial

test_parcel_moving_equidistant_SOURCES = test_utils.f90 test_parcel_moving_equidistant.f90
test_parcel_moving_equidistant_LDADD = libcombi.la

test_parcel_moving_random_SOURCES = test_utils.f90 test_parcel_moving_random.f90
test_parcel_moving_random_LDADD = libcombi.la

test_parcel_split_random_SOURCES = test_utils.f90 test_parcel_split_random.f90
test_parcel_split_random_LDADD = libcombi.la

test_parcel_merge_random_SOURCES = test_utils.f90 test_parcel_merge_random.f90
test_parcel_merge_random_LDADD = libcombi.la

test_parcel_split_merge_SOURCES = test_utils.f90 test_parcel_split_merge.f90
test_parcel_split_merge_LDADD = libcombi.la

test_merging_parcels_SOURCES = test_utils.f90 test_merging_parcels.f90
test_merging_parcels_LDADD = libcombi.la

test_merging_parcels_serial_SOURCES = test_utils.f90			\
				      parcel_nearest_serial.f90		\
				      parcel_merge_serial.f90		\
				      test_merging_parcels_serial.f90
test_merging_parcels_serial_LDADD = libcombi.la

nobase_mpitests_DATA = test_nearest_algorithm.py

endif
