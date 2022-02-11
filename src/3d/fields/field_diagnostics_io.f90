module field_diagnostics_io
    use field_diagnostics_hdf5
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

            call start_timer(field_stat_io_timer)

            call create_h5_field_stats_file(basename, overwrite, l_restart)

#ifdef ENABLE_NETCDF
            call create_netcdf_field_stats_file(basename, overwrite, l_restart)
#endif
            call stop_timer(field_stat_io_timer)
        end subroutine create_field_stats_file

        subroutine write_field_stats_step(t, dt)
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            call start_timer(field_stat_io_timer)

            call calculate_field_diagnostics

            call write_h5_field_stats_step(t, dt)

#ifdef ENABLE_NETCDF
            call write_h5_field_stats_step(t, dt)
#endif
            call stop_timer(field_stat_io_timer)

        end subroutine write_field_stats_step

end module field_diagnostics_io
