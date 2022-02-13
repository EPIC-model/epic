module parcel_io
    use parcel_netcdf
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: parcel_io_timer

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        ! @param[in] l_restart if in restart mode
        subroutine create_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

            call create_netcdf_parcel_file(basename,    &
                                           overwrite,   &
                                           l_restart)
        end subroutine create_parcel_file

        ! Write parcels
        ! @param[in] t is the time
        subroutine write_parcels(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_io_timer)

            call write_netcdf_parcels(t)
            call stop_timer(parcel_io_timer)
        end subroutine write_parcels

        ! Read parcels
        subroutine read_parcels(fname)
            character(*), intent(in) :: fname

            call start_timer(parcel_io_timer)

            call read_netcdf_parcels(fname)

            call stop_timer(parcel_io_timer)
        end subroutine read_parcels

end module parcel_io
