AM_FCFLAGS =				\
	-I $(top_builddir)/src/3d/ 	\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =						\
	$(top_builddir)/src/3d/parcels/*.o 		\
	$(top_builddir)/src/3d/utils/*.o		\
	$(top_builddir)/src/3d/fields/*.o		\
	$(top_builddir)/src/3d/inversion/*.o

# 21 November 2021
# gnu.org/software/automake/manual/html_node/Libtool-Convenience-Libraries.html
noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/mpi/libepic_mpi.la		\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/fft/libepic_mpi_fft.la	\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

unittestsdir = $(prefix)/unit-tests

# 	test_ellipsoid_split
# 	test_ellipsoid_reflection
# 	test_free_slip_3d
unittests_PROGRAMS = 					\
	test_jacobi_1					\
	test_jacobi_2					\
	test_jacobi_3					\
	test_jacobi_4					\
	test_scherzinger_1				\
	test_scherzinger_2				\
	test_scherzinger_3				\
	test_scherzinger_4				\
	test_mpi_parcel_init_3d				\
	test_mpi_trilinear				\
	test_mpi_laplace_correction_3d			\
	test_mpi_gradient_correction_3d			\
	test_mpi_parcel_correction_3d			\
	test_mpi_grid2par				\
	test_mpi_field_halo_fill			\
	test_mpi_field_halo_fill_vector			\
	test_mpi_field_interior_accumulate		\
	test_mpi_field_interior_accumulate_vector	\
	test_mpi_field_diagnostics			\
	test_mpi_parcel_diagnostics			\
	test_mpi_parcel_write				\
	test_mpi_parcel_delete				\
	test_mpi_parcel_pack				\
	test_mpi_parcel_unpack				\
	test_mpi_parcel_read				\
	test_mpi_parcel_read_rejection			\
	test_mpi_parcel_communicate			\
	test_mpi_parcel_split				\
	test_mpi_fft					\
	test_mpi_dst					\
	test_mpi_reverse_x				\
	test_mpi_reverse_y				\
	test_mpi_diffx					\
	test_mpi_diffy					\
	test_mpi_diffz0					\
	test_mpi_diffz1					\
	test_mpi_diffz2					\
	test_mpi_diffz3					\
 	test_mpi_diffz4					\
 	test_mpi_diffz5					\
 	test_mpi_fft0					\
 	test_mpi_fft1					\
 	test_mpi_vor2vel				\
 	test_mpi_vel2vgrad				\
	test_mpi_vtend					\
	test_mpi_nearest_1				\
	test_mpi_nearest_2				\
	test_mpi_nearest_3				\
	test_mpi_nearest_4				\
	test_mpi_nearest_5				\
	test_mpi_nearest_6				\
	test_mpi_nearest_7				\
	test_mpi_nearest_8				\
	test_mpi_nearest_9				\
	test_mpi_nearest_10				\
	test_mpi_nearest_11				\
	test_mpi_nearest_12				\
	test_mpi_nearest_13				\
	test_mpi_nearest_14				\
	test_mpi_nearest_15				\
	test_mpi_nearest_16				\
	test_mpi_nearest_17				\
	test_mpi_nearest_18				\
	test_mpi_nearest_19				\
	test_mpi_nearest_20
# 	test_mpi_lapinv1

dataroot_DATA =		\
	numpy_eigh.py 	\
	V.asc	   	\
	D.asc	   	\
	B.asc

test_jacobi_1_SOURCES = test_jacobi_1.f90
test_jacobi_1_LDADD = libcombi.la

test_jacobi_2_SOURCES = test_jacobi_2.f90
test_jacobi_2_LDADD = libcombi.la

test_jacobi_3_SOURCES = test_jacobi_3.f90
test_jacobi_3_LDADD = libcombi.la

test_jacobi_4_SOURCES = test_jacobi_4.f90
test_jacobi_4_LDADD = libcombi.la

test_scherzinger_1_SOURCES = test_scherzinger_1.f90
test_scherzinger_1_LDADD = libcombi.la

test_scherzinger_2_SOURCES = test_scherzinger_2.f90
test_scherzinger_2_LDADD = libcombi.la

test_scherzinger_3_SOURCES = test_scherzinger_3.f90
test_scherzinger_3_LDADD = libcombi.la

test_scherzinger_4_SOURCES = test_scherzinger_4.f90
test_scherzinger_4_LDADD = libcombi.la

# test_ellipsoid_split_SOURCES = test_ellipsoid_split.f90
# test_ellipsoid_split_LDADD = libcombi.la

test_mpi_parcel_init_3d_SOURCES = test_mpi_parcel_init_3d.f90
test_mpi_parcel_init_3d_LDADD = libcombi.la

test_mpi_trilinear_SOURCES = test_mpi_trilinear.f90
test_mpi_trilinear_LDADD = libcombi.la

test_mpi_laplace_correction_3d_SOURCES = test_mpi_laplace_correction_3d.f90
test_mpi_laplace_correction_3d_LDADD = libcombi.la

test_mpi_gradient_correction_3d_SOURCES = test_mpi_gradient_correction_3d.f90
test_mpi_gradient_correction_3d_LDADD = libcombi.la

test_mpi_parcel_correction_3d_SOURCES = test_mpi_parcel_correction_3d.f90
test_mpi_parcel_correction_3d_LDADD = libcombi.la

# test_ellipsoid_reflection_SOURCES = test_ellipsoid_reflection.f90
# test_ellipsoid_reflection_LDADD = libcombi.la

test_mpi_grid2par_SOURCES = test_mpi_grid2par.f90
test_mpi_grid2par_LDADD = libcombi.la

# test_par2grid_SOURCES = test_par2grid.f90
# test_par2grid_LDADD = libcombi.la

# test_free_slip_3d_SOURCES = test_free_slip_3d.f90
# test_free_slip_3d_LDADD = libcombi.la

test_mpi_field_halo_fill_SOURCES = test_mpi_field_halo_fill.f90
test_mpi_field_halo_fill_LDADD = libcombi.la

test_mpi_field_halo_fill_vector_SOURCES = test_mpi_field_halo_fill_vector.f90
test_mpi_field_halo_fill_vector_LDADD = libcombi.la

test_mpi_field_interior_accumulate_SOURCES = test_mpi_field_interior_accumulate.f90
test_mpi_field_interior_accumulate_LDADD = libcombi.la

test_mpi_field_interior_accumulate_vector_SOURCES = test_mpi_field_interior_accumulate_vector.f90
test_mpi_field_interior_accumulate_vector_LDADD = libcombi.la

test_mpi_field_diagnostics_SOURCES = test_mpi_field_diagnostics.f90
test_mpi_field_diagnostics_LDADD = libcombi.la

test_mpi_parcel_diagnostics_SOURCES = test_mpi_parcel_diagnostics.f90
test_mpi_parcel_diagnostics_LDADD = libcombi.la

test_mpi_parcel_write_SOURCES = test_mpi_parcel_write.f90
test_mpi_parcel_write_LDADD = libcombi.la

test_mpi_parcel_delete_SOURCES = test_mpi_parcel_delete.f90
test_mpi_parcel_delete_LDADD = libcombi.la

test_mpi_parcel_pack_SOURCES = test_mpi_parcel_pack.f90
test_mpi_parcel_pack_LDADD = libcombi.la

test_mpi_parcel_unpack_SOURCES = test_mpi_parcel_unpack.f90
test_mpi_parcel_unpack_LDADD = libcombi.la

test_mpi_parcel_read_SOURCES = test_mpi_parcel_read.f90
test_mpi_parcel_read_LDADD = libcombi.la

test_mpi_parcel_read_rejection_SOURCES = test_mpi_parcel_read_rejection.f90
test_mpi_parcel_read_rejection_LDADD = libcombi.la

test_mpi_parcel_communicate_SOURCES = test_mpi_parcel_communicate.f90
test_mpi_parcel_communicate_LDADD = libcombi.la

test_mpi_parcel_split_SOURCES = test_mpi_parcel_split.f90
test_mpi_parcel_split_LDADD = libcombi.la

test_mpi_fft_SOURCES = test_mpi_fft.f90
test_mpi_fft_LDADD = libcombi.la

test_mpi_dst_SOURCES = test_mpi_dst.f90
test_mpi_dst_LDADD = libcombi.la

test_mpi_reverse_x_SOURCES = test_mpi_reverse_x.f90
test_mpi_reverse_x_LDADD = libcombi.la

test_mpi_reverse_y_SOURCES = test_mpi_reverse_y.f90
test_mpi_reverse_y_LDADD = libcombi.la

test_mpi_diffx_SOURCES = test_mpi_diffx.f90
test_mpi_diffx_LDADD = libcombi.la

test_mpi_diffy_SOURCES = test_mpi_diffy.f90
test_mpi_diffy_LDADD = libcombi.la

test_mpi_diffz0_SOURCES = test_mpi_diffz0.f90
test_mpi_diffz0_LDADD = libcombi.la

test_mpi_diffz1_SOURCES = test_mpi_diffz1.f90
test_mpi_diffz1_LDADD = libcombi.la

test_mpi_diffz2_SOURCES = test_mpi_diffz2.f90
test_mpi_diffz2_LDADD = libcombi.la

test_mpi_diffz3_SOURCES = test_mpi_diffz3.f90
test_mpi_diffz3_LDADD = libcombi.la

test_mpi_diffz4_SOURCES = test_mpi_diffz4.f90
test_mpi_diffz4_LDADD = libcombi.la

test_mpi_diffz5_SOURCES = test_mpi_diffz5.f90
test_mpi_diffz5_LDADD = libcombi.la

test_mpi_fft0_SOURCES = test_mpi_fft0.f90
test_mpi_fft0_LDADD = libcombi.la

test_mpi_fft1_SOURCES = test_mpi_fft1.f90
test_mpi_fft1_LDADD = libcombi.la

test_mpi_vor2vel_SOURCES = test_mpi_vor2vel.f90
test_mpi_vor2vel_LDADD = libcombi.la

test_mpi_vel2vgrad_SOURCES = test_mpi_vel2vgrad.f90
test_mpi_vel2vgrad_LDADD = libcombi.la

test_mpi_vtend_SOURCES = test_mpi_vtend.f90
test_mpi_vtend_LDADD = libcombi.la

test_mpi_nearest_1_SOURCES = test_mpi_nearest_1.f90
test_mpi_nearest_1_LDADD = libcombi.la

test_mpi_nearest_2_SOURCES = test_mpi_nearest_2.f90
test_mpi_nearest_2_LDADD = libcombi.la

test_mpi_nearest_3_SOURCES = test_mpi_nearest_3.f90
test_mpi_nearest_3_LDADD = libcombi.la

test_mpi_nearest_4_SOURCES = test_mpi_nearest_4.f90
test_mpi_nearest_4_LDADD = libcombi.la

test_mpi_nearest_5_SOURCES = test_mpi_nearest_5.f90
test_mpi_nearest_5_LDADD = libcombi.la

test_mpi_nearest_6_SOURCES = test_mpi_nearest_6.f90
test_mpi_nearest_6_LDADD = libcombi.la

test_mpi_nearest_7_SOURCES = test_mpi_nearest_7.f90
test_mpi_nearest_7_LDADD = libcombi.la

test_mpi_nearest_8_SOURCES = test_mpi_nearest_8.f90
test_mpi_nearest_8_LDADD = libcombi.la

test_mpi_nearest_9_SOURCES = test_mpi_nearest_9.f90
test_mpi_nearest_9_LDADD = libcombi.la

test_mpi_nearest_10_SOURCES = test_mpi_nearest_10.f90
test_mpi_nearest_10_LDADD = libcombi.la

test_mpi_nearest_11_SOURCES = test_mpi_nearest_11.f90
test_mpi_nearest_11_LDADD = libcombi.la

test_mpi_nearest_12_SOURCES = test_mpi_nearest_12.f90
test_mpi_nearest_12_LDADD = libcombi.la

test_mpi_nearest_13_SOURCES = test_mpi_nearest_13.f90
test_mpi_nearest_13_LDADD = libcombi.la

test_mpi_nearest_14_SOURCES = test_mpi_nearest_14.f90
test_mpi_nearest_14_LDADD = libcombi.la

test_mpi_nearest_15_SOURCES = test_mpi_nearest_15.f90
test_mpi_nearest_15_LDADD = libcombi.la

test_mpi_nearest_16_SOURCES = test_mpi_nearest_16.f90
test_mpi_nearest_16_LDADD = libcombi.la

test_mpi_nearest_17_SOURCES = test_mpi_nearest_17.f90
test_mpi_nearest_17_LDADD = libcombi.la

test_mpi_nearest_18_SOURCES = test_mpi_nearest_18.f90
test_mpi_nearest_18_LDADD = libcombi.la

test_mpi_nearest_19_SOURCES = test_mpi_nearest_19.f90
test_mpi_nearest_19_LDADD = libcombi.la

test_mpi_nearest_20_SOURCES = test_mpi_nearest_20.f90
test_mpi_nearest_20_LDADD = libcombi.la

# test_mpi_lapinv1_SOURCES = test_mpi_lapinv1.f90
# test_mpi_lapinv1_LDADD = libcombi.la
