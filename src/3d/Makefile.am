if ENABLE_3D
AM_FCFLAGS =				\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/netcdf

bin_PROGRAMS = epic3d
epic3d_SOURCES = 				\
	utils/dimensions.f90			\
	utils/options.f90			\
	utils/parameters.f90			\
	utils/interpl.f90			\
	fields/field_mpi.f90			\
	fields/fields.f90			\
	fields/field_ops.f90			\
	fields/field_diagnostics.f90		\
	fields/field_netcdf.f90			\
	fields/field_diagnostics_netcdf.f90	\
	parcels/parcel_ops.f90			\
	parcels/parcel_container.f90 		\
	parcels/parcel_ellipsoid.f90		\
	parcels/parcel_ellipse.f90		\
	parcels/parcel_ellipse_init.f90		\
	parcels/parcel_ellipse_split.f90	\
	parcels/parcels.f90			\
	parcels/parcel_bc.f90                   \
	parcels/parcel_mpi.f90			\
	parcels/parcel_nearest.f90		\
	parcels/parcel_merge.f90		\
	parcels/parcel_split.f90		\
	parcels/parcel_diagnostics.f90		\
	parcels/parcel_ellipse_interpl.f90	\
	parcels/parcel_interpl.f90 		\
	parcels/parcel_init.f90			\
	inversion/inversion_utils.f90		\
	inversion/inversion.f90			\
	parcels/parcel_correction.f90		\
	parcels/parcel_netcdf.f90		\
	parcels/parcel_diagnostics_netcdf.f90	\
	boundary_layer/bndry_fluxes.f90		\
	utils/utils.f90				\
	stepper/rk_utils.f90			\
	parcels/parcel_damping.f90	        \
	stepper/ls_rk.f90			\
	epic3d.f90

epic3d_LDADD = 						\
	$(top_builddir)/src/utils/libepic_utils.la 	\
	$(top_builddir)/src/fft/libepic_mpi_fft.la	\
	$(top_builddir)/src/mpi/libepic_mpi.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod *.smod
endif
