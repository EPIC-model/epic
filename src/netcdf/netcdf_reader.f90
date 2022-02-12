! =============================================================================
!                               NetCDF reader
!   This module reads in datatsets. For further information about the
!   NetCDF Fortran API please visit
!   https://www.unidata.ucar.edu/software/netcdf/docs-fortran/index.html
! =============================================================================
module netcdf_reader
    use netcdf_utils
    implicit none

    interface read_netcdf_global_attribute
        module procedure :: read_netcdf_global_attrib_integer
        module procedure :: read_netcdf_global_attrib_double
    end interface read_netcdf_global_attribute

    private :: read_netcdf_global_attrib_integer,   &
               read_netcdf_global_attrib_double

    contains

        ! @returns -1 if NetCDF file does not contain the attribute 'n_parcels'
        subroutine get_num_parcels(ncid, n_parcels)
            integer, intent(in)  :: ncid
            integer, intent(out) :: n_parcels

            if (.not. has_attribute(ncid, 'n_parcels')) then
                n_parcels = -1
                return
            endif

            call read_netcdf_global_attrib_integer(ncid, 'n_parcels', n_parcels)
        end subroutine get_num_parcels

        ! @returns -1 if NetCDF file does not contain the attribute 'n_steps'
        subroutine get_num_steps(ncid, n_steps)
            integer, intent(in)  :: ncid
            integer, intent(out) :: n_steps

            if (.not. has_attribute(ncid, 'n_steps')) then
                n_steps = -1
                return
            endif

            call read_netcdf_global_attrib_integer(ncid, 'n_steps', n_steps)
        end subroutine get_num_steps

        ! @returns -1 if NetCDF file does not contain the attribute 't'
        subroutine get_time(ncid, t)
            integer,          intent(in)  :: ncid
            double precision, intent(out) :: t

            if (has_attribute(ncid, 't')) then
                call read_netcdf_global_attrib_double(ncid, 't', t)
                return
            endif

!             if (has_dataset(ncid, 't')) then
!                 call read_netcdf_dataset(ncid, 't', t)
!                 return
!             endif

            t = -1.0d0
        end subroutine get_time

        ! 11 Jan 2022
        ! https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/FORTRAN/H5T/h5ex_t_stringCatt_F03.f90
        subroutine get_file_type(ncid, file_type)
            integer, intent(in)  :: ncid
            character(:), allocatable, intent(out) :: file_type

            if (.not. has_attribute(ncid, 'file_type')) then
                print *, 'Not a proper EPIC NetCDF file.'
                stop
            endif

        end subroutine get_file_type

        function has_attribute(ncid, name) result(link_exists)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            logical                  :: link_exists
            link_exists = .false.

            ncerr = nf90_inquire_attribute(ncid, NF90_GLOBAL, name)

            if (ncerr /= nf90_noerr) then
                link_exists = .true.
            endif
            ncerr = 0
        end function has_attribute

        function has_dataset(ncid, name) result(link_exists)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            integer                  :: varid
            logical                  :: link_exists
            link_exists = .false.

            ncerr = nf90_inq_varid(ncid, name, varid)
            ncerr = nf90_inquire_variable(ncid, varid)
            if (ncerr /= nf90_noerr) then
                link_exists = .true.
            endif
            ncerr = 0
        end function has_dataset

        subroutine read_netcdf_dataset(ncid, name, buffer)
            integer,          intent(in)  :: ncid
            character(*),     intent(in)  :: name
            double precision, intent(out) :: buffer(:)
            integer                       :: varid

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading dataset id failed.")

!             ncerr = nf90_get_var(ncid, varid, values, start, count

        end subroutine read_netcdf_dataset

        subroutine read_netcdf_global_attrib_integer(ncid, name, val)
            integer,       intent(in)     :: ncid
            character(*),  intent(in)     :: name
            integer,       intent(out)    :: val

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, name, val)
            call check_netcdf_error("Reading attribute '" // name // "' failed.")
        end subroutine read_netcdf_global_attrib_integer

        subroutine read_netcdf_global_attrib_double(ncid, name, val)
            integer,          intent(in)     :: ncid
            character(*),     intent(in)     :: name
            double precision, intent(out)    :: val

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, name, val)
            call check_netcdf_error("Reading attribute '" // name // "' failed.")

        end subroutine read_netcdf_global_attrib_double

        subroutine read_netcdf_domain(ncfname, origin, extent, ncells)
            character(*), intent(in)      :: ncfname
            integer                       :: ncid
            double precision, intent(out) :: extent(:), origin(:)
            integer,          intent(out) :: ncells(:)

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
            call get_netcdf_box(ncid, origin, extent, ncells)
            call close_netcdf_file(ncid)

        end subroutine read_netcdf_domain

        subroutine get_netcdf_box(ncid, origin, extent, ncells)
            integer,          intent(in)     :: ncid
            double precision, intent(out)    :: extent(:), origin(:)
            integer,          intent(out)    :: ncells(:)

            if ((size(ncells) > 3) .or. (size(extent) > 3) .or. (size(extent) > 3)) then
                print *, "Cannot read more than 3 dimensions!"
                stop
            endif

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, "ncells", ncells)
            call check_netcdf_error("Reading attribute 'ncells' failed.")

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, "extent", extent)
            call check_netcdf_error("Reading attribute 'extent' failed.")

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, "origin", origin)
            call check_netcdf_error("Reading attribute 'origin' failed.")

        end subroutine get_netcdf_box

end module netcdf_reader
