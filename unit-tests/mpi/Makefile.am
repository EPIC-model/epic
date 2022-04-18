AM_FCFLAGS =				\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/3d/	\
	-I $(top_builddir)/src/mpi

noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/mpi/libepic_mpi.la

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_mpi_init				\
	test_mpi_layout

test_mpi_init_SOURCES = test_mpi_init.f90
test_mpi_init_LDADD = libcombi.la

test_mpi_layout_SOURCES = test_mpi_layout.f90
test_mpi_layout_LDADD = libcombi.la
