bin_PROGRAMS = epic
epic_SOURCES = 					\
	utils/constants.f90 			\
	utils/jacobi.f90			\
	utils/timer.f90				\
	utils/phys_constants.f90		\
	utils/phys_parameters.f90		\
	utils/stafft.f90			\
	utils/sta2dfft.f90			\
	utils/deriv1d.f90			\
	utils/options.f90			\
	utils/merge_sort.f90			\
	utils/parser.f90 			\
	utils/parameters.f90 			\
	hdf5/h5_utils.f90			\
	hdf5/h5_writer.f90			\
	hdf5/h5_reader.f90			\
	fields/fields.f90			\
	fields/field_diagnostics.f90		\
	fields/field_hdf5.f90			\
	parcels/parcel_ellipse.f90		\
	parcels/parcel_ellipsoid.f90		\
	parcels/parcel_container.f90 		\
	parcels/parcel_bc.f90			\
	parcels/parcel_nearest.f90		\
	parcels/merge_hdf5.f90			\
	parcels/parcel_merge.f90		\
	parcels/parcel_split.f90		\
	parcels/parcel_diagnostics.f90		\
	parcels/parcel_interpl.f90 		\
	parcels/parcel_init.f90			\
	parcels/parcel_correction.f90		\
	parcels/parcel_hdf5.f90			\
	tri_inversion.f90			\
	utils/utils.f90				\
	stepper/rk4_utils.f90			\
	stepper/ls_rk4.f90			\
	epic.f90

clean-local:
	rm -f *.mod
