module parcel_io
#ifdef ENABLE_NETCDF
    use parcel_netcdf
#else
    use parcel_hdf5
#endif
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

#ifdef ENABLE_NETCDF
            call create_netcdf_parcel_file(basename,    &
                                           overwrite,   &
                                           l_restart)
#else
            call create_h5_parcel_file(basename,        &
                                       overwrite,       &
                                       l_restart)
#endif
        end subroutine create_parcel_file

        ! Write parcels
        ! @param[in] t is the time
        subroutine write_parcels(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_io_timer)

#ifdef ENABLE_NETCDF
            call write_netcdf_parcels(t)
#else
            call write_h5_parcels(t)
#endif
            call stop_timer(parcel_io_timer)
        end subroutine write_parcels

        ! Read parcels
        subroutine read_parcels(fname)
            character(*), intent(in) :: fname

            call start_timer(parcel_io_timer)

#ifdef ENABLE_NETCDF
            call read_netcdf_parcels(fname)
#else
            call read_h5_parcels(fname)
#endif
            call stop_timer(parcel_io_timer)
        end subroutine read_parcels

end module parcel_io
