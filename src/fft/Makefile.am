lib_LTLIBRARIES = libepic_fft.la

if ENABLE_3D
lib_LTLIBRARIES += libepic_mpi_fft.la
AM_FCFLAGS = -I $(top_builddir)/src/mpi
endif

libepic_fft_la_SOURCES = 	\
	../utils/constants.f90	\
	stafft.f90		\
	sta2dfft.f90

if ENABLE_3D
libepic_mpi_fft_la_SOURCES = 		\
	../utils/constants.f90		\
	../3d/utils/dimensions.f90	\
	../utils/datatypes.f90		\
	stafft.f90			\
	sta2dfft.f90			\
	pencilfft.f90

libepic_mpi_fft_la_LIBADD = $(top_builddir)/src/mpi/libepic_mpi.la
endif

