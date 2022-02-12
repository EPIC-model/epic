module parcel_diagnostics_io
#ifdef ENABLE_HDF5
    use parcel_diagnostics_hdf5
#endif
#ifdef ENABLE_NETCDF
    use parcel_diagnostics_netcdf
#endif
    use parcel_diagnostics
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: parcel_stat_io_timer

    contains

        subroutine create_parcel_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

#ifdef ENABLE_HDF5
            call create_h5_parcel_stats_file(basename, overwrite, l_restart)
#endif

#ifdef ENABLE_NETCDF
            call create_netcdf_parcel_stats_file(basename, overwrite, l_restart)
#endif
        end subroutine create_parcel_stats_file

        subroutine write_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_stat_io_timer)

#ifdef ENABLE_HDF5
            call write_h5_parcel_stats(t)
#endif

#ifdef ENABLE_NETCDF
            call write_netcdf_parcel_stats(t)
#endif
            call stop_timer(parcel_stat_io_timer)

        end subroutine write_parcel_stats

end module parcel_diagnostics_io
