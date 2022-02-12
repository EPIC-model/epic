module field_diagnostics_io
#ifdef ENABLE_NETCDF
    use field_diagnostics_netcdf
#else
    use field_diagnostics_hdf5
#endif
    use field_diagnostics
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_stat_io_timer

    public :: create_field_stats_file,  &
              write_field_stats

    contains

        subroutine create_field_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

#ifdef ENABLE_NETCDF
            call create_netcdf_field_stats_file(basename, overwrite, l_restart)
#else
            call create_h5_field_stats_file(basename, overwrite, l_restart)
#endif
        end subroutine create_field_stats_file

        subroutine write_field_stats(t)
            double precision, intent(in)    :: t

            call start_timer(field_stat_io_timer)

            call calculate_field_diagnostics

#ifdef ENABLE_NETCDF
            call write_netcdf_field_stats(t)
#else
            call write_h5_field_stats(t)
#endif
            call stop_timer(field_stat_io_timer)

        end subroutine write_field_stats

end module field_diagnostics_io
