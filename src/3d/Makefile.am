if ENABLE_3D
AM_FCFLAGS =				\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/netcdf

bin_PROGRAMS = epic3d
epic3d_SOURCES = 				\
	config.f90				\
	utils/deriv1d.f90			\
	utils/options.f90			\
	utils/parameters.f90			\
	utils/phys_parameters.f90		\
	fields/fields.f90			\
	fields/field_diagnostics.f90		\
	fields/field_netcdf.f90			\
	fields/field_diagnostics_netcdf.f90	\
	fields/field_io.f90			\
	fields/field_diagnostics_io.f90		\
	parcels/parcel_ellipsoid.f90		\
	parcels/parcel_container.f90 		\
	parcels/parcel_bc.f90			\
	parcels/parcel_nearest.f90		\
	parcels/parcel_merge.f90		\
	parcels/parcel_split.f90		\
	parcels/parcel_diagnostics.f90		\
	parcels/parcel_interpl.f90 		\
	parcels/parcel_init.f90			\
	inversion/inversion_utils.f90		\
	inversion/inversion.f90			\
	parcels/parcel_correction.f90		\
	parcels/parcel_netcdf.f90		\
	parcels/parcel_diagnostics_netcdf.f90	\
	parcels/parcel_io.f90			\
	parcels/parcel_diagnostics_io.f90	\
	utils/utils.f90				\
	stepper/rk4_utils.f90			\
	stepper/ls_rk4.f90			\
	epic3d.f90

epic3d_LDADD = 						\
	$(top_builddir)/src/utils/libepic_utils.la 	\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod
endif
