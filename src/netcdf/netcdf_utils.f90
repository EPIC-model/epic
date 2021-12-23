! =============================================================================
! References:
! https://github.com/Unidata/netcdf-fortran/tree/master/examples (23.12.2021)
! https://docs.unidata.ucar.edu/netcdf-fortran/current/f90_The-NetCDF-Fortran-90-Interface-Guide.html
! =============================================================================
module netcdf_utils
    use netcdf
    implicit none

    ! netCDF error if non-zero
    integer :: ncerr = 0

    contains

        subroutine create_netcdf_file(ncfname, overwrite, ncid)
            character(*), intent(in)  :: ncfname
            logical,      intent(in)  :: overwrite
            integer,      intent(out) :: ncid
            logical                     :: exists = .true.

            ! check whether file exists
            inquire(file=ncfname, exist=exists)

            if (exists .and. overwrite) then
                call delete_netcdf_file(ncfname)
            else if (exists) then
                print *, "File '" // trim(ncfname) // "' already exists. Exiting."
                stop
            endif

            ! nf90_noclobber -- do not want to clobber (overwrite) an existing dataset
            ! nf90_clobber -- overwrite this file, if it already exists.
            ncerr = nf90_create(path = ncfname,        &
                                cmode = nf90_clobber,  &
                                ncid = ncid)

            call check_netcdf_error("Failed to create netcdf file'" // trim(ncfname) // "'.")

        end subroutine create_netcdf_file

        subroutine delete_netcdf_file(ncfname)
            character(*), intent(in) :: ncfname
            integer                  :: stat

            ! 15 June 2021
            ! https://stackoverflow.com/questions/18668832/how-delete-file-from-fortran-code
            open(unit=1234, iostat=stat, file=ncfname, status='old')
            if (stat == 0) then
                close(1234, status='delete')
            endif
        end subroutine delete_netcdf_file


        subroutine open_netcdf_file(ncfname, access_flag, ncid)
            character(*), intent(in)  :: ncfname
            integer,      intent(in)  :: access_flag ! NF90_WRITE or NF90_NOWRITE
            integer,      intent(out) :: ncid

            ncerr = nf90_open(path = ncfname,       &
                              mode = access_flag,   &
                              ncid = ncid)

            call check_netcdf_error("Opening the netcdf file failed.")

        end subroutine open_netcdf_file


        subroutine close_netcdf_file(ncid)
            integer, intent(in) :: ncid

            ncerr = nf90_close(ncid)

            call check_netcdf_error("Closing the netcdf file failed.")
        end subroutine close_netcdf_file
!
        subroutine check_netcdf_error(msg)
            character(*), intent(in) :: msg
#ifndef NDEBUG
            if (ncerr /= nf90_noerr) then
                print *, msgq
                print *, trim(nf90_strerror(ncerr))
                stop
            endif
#endif
        end subroutine check_netcdf_error
end module netcdf_utils
