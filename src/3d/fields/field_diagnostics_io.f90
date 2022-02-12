module field_diagnostics_io
#ifdef ENABLE_HDF5
    use field_diagnostics_hdf5
#endif
#ifdef ENABLE_NETCDF
    use field_diagnostics_netcdf
#endif
    use field_diagnostics
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_stat_io_timer

    public :: create_field_stats_file,  &
              write_field_stats_step

    contains

        subroutine create_field_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

#ifdef ENABLE_HDF5
            call create_h5_field_stats_file(basename, overwrite, l_restart)
#endif

#ifdef ENABLE_NETCDF
            call create_netcdf_field_stats_file(basename, overwrite, l_restart)
#endif
        end subroutine create_field_stats_file

        subroutine write_field_stats_step(t, dt)
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            call start_timer(field_stat_io_timer)

            call calculate_field_diagnostics

#ifdef ENABLE_HDF5
            call write_h5_field_stats_step(t, dt)
#endif

#ifdef ENABLE_NETCDF
            call write_netcdf_field_stats_step(t)
#endif
            call stop_timer(field_stat_io_timer)

        end subroutine write_field_stats_step

end module field_diagnostics_io
