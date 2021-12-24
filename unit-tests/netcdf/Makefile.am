if ENABLE_UNIT_TESTS
if ENABLE_NETCDF

AM_FCFLAGS =				\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/utils    \
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =                                            \
        $(top_builddir)/src/3d/utils/*.o

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_netcdf_dataset_2d			\
	test_netcdf_dataset_3d			\
	test_netcdf_time_series_1		\
	test_netcdf_time_series_2

test_netcdf_dataset_2d_SOURCES = test_netcdf_dataset_2d.f90
test_netcdf_dataset_2d_LDADD = libcombi.la

test_netcdf_dataset_3d_SOURCES = test_netcdf_dataset_3d.f90
test_netcdf_dataset_3d_LDADD = libcombi.la

test_netcdf_time_series_1_SOURCES = test_netcdf_time_series_1.f90
test_netcdf_time_series_1_LDADD = libcombi.la

test_netcdf_time_series_2_SOURCES = test_netcdf_time_series_2.f90
test_netcdf_time_series_2_LDADD = libcombi.la

endif
endif
