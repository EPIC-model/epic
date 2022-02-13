module parcel_diagnostics_io
    use parcel_diagnostics_netcdf
    use parcel_diagnostics
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: parcel_stat_io_timer

    contains

        subroutine create_parcel_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

            call create_netcdf_parcel_stats_file(basename, overwrite, l_restart)
        end subroutine create_parcel_stats_file

        subroutine write_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_stat_io_timer)

            call write_netcdf_parcel_stats(t)

            call stop_timer(parcel_stat_io_timer)
        end subroutine write_parcel_stats

end module parcel_diagnostics_io
