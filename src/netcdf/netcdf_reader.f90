! =============================================================================
!                               NetCDF reader
!   This module reads in datatsets. For further information about the
!   NetCDF Fortran API please visit
!   https://www.unidata.ucar.edu/software/netcdf/docs-fortran/index.html
! =============================================================================
module netcdf_reader
    use netcdf_utils
    implicit none

    interface read_netcdf_dataset
        module procedure :: read_netcdf_dataset_1d
        module procedure :: read_netcdf_dataset_2d
        module procedure :: read_netcdf_dataset_3d
        module procedure :: read_netcdf_dataset_1d_integer
    end interface read_netcdf_dataset

    interface read_netcdf_attribute
        module procedure :: read_netcdf_attrib_integer
        module procedure :: read_netcdf_attrib_double
        module procedure :: read_netcdf_attrib_character
        module procedure :: read_netcdf_attrib_logical
    end interface read_netcdf_attribute

    interface read_netcdf_attribute_default
        module procedure :: read_netcdf_attrib_default_integer
        module procedure :: read_netcdf_attrib_default_double
        module procedure :: read_netcdf_attrib_default_logical
    end interface read_netcdf_attribute_default

    private :: read_netcdf_attrib_integer,          &
               read_netcdf_attrib_double,           &
               read_netcdf_attrib_character,        &
               read_netcdf_attrib_logical,          &
               read_netcdf_attrib_default_integer,  &
               read_netcdf_attrib_default_double,   &
               read_netcdf_attrib_default_logical,  &
               read_netcdf_dataset_1d,              &
               read_netcdf_dataset_2d,              &
               read_netcdf_dataset_3d

    contains

        subroutine get_var_id(ncid, name, varid)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            integer,      intent(out) :: varid

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading the ID of '" // name // "' failed.")
        end subroutine get_var_id

        subroutine get_dim_id(ncid, name, dimid)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            integer,      intent(out) :: dimid

            ncerr = nf90_inq_dimid(ncid, name, dimid)
            call check_netcdf_error("Reading the ID of '" // name // "' failed.")
        end subroutine get_dim_id

        subroutine get_dimension_size(ncid, name, num)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            integer, intent(out)     :: num
            integer                  :: dimid

            ncerr = nf90_inq_dimid(ncid, name, dimid)
            call check_netcdf_error("Reading " // name // " dimension id failed.")
            ncerr = nf90_inquire_dimension(ncid, dimid, len=num)
            call check_netcdf_error("Reading " // name // " failed.")
        end subroutine get_dimension_size

        subroutine get_num_parcels(ncid, n_parcels)
            integer, intent(in)  :: ncid
            integer, intent(out) :: n_parcels

            call get_dimension_size(ncid, 'n_parcels', n_parcels)
        end subroutine get_num_parcels

        subroutine get_num_steps(ncid, n_steps)
            integer, intent(in)  :: ncid
            integer, intent(out) :: n_steps

            call get_dimension_size(ncid, 't', n_steps)
        end subroutine get_num_steps

        subroutine get_time(ncid, t)
            integer,          intent(in)  :: ncid
            double precision, intent(out) :: t
            integer                       :: n_steps, varid, start(1), cnt(1)
            double precision              :: values(1)

            call get_num_steps(ncid, n_steps)

            if (has_dataset(ncid, 't')) then
                start(1) = n_steps
                cnt(1) = 1
                ncerr = nf90_inq_varid(ncid, 't', varid)
                call check_netcdf_error("Reading time id failed.")
                ncerr = nf90_get_var(ncid, varid, values, start=start, count=cnt)
                t = values(1)
                call check_netcdf_error("Reading time failed.")
                return
            else
                print *, "Error: No time dataset found."
                stop
            endif
        end subroutine get_time

        subroutine get_file_type(ncid, file_type)
            integer,      intent(in)  :: ncid
            character(*), intent(out) :: file_type

            if (.not. has_attribute(ncid, 'file_type')) then
                print *, 'Not a proper EPIC NetCDF file.'
                stop
            endif

            ncerr = nf90_get_att(ncid=ncid, varid=NF90_GLOBAL, &
                                 name='file_type', values=file_type)
            call check_netcdf_error("Reading file type failed.")

        end subroutine get_file_type

        function has_attribute(ncid, name) result(link_exists)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            logical                  :: link_exists

            ncerr = nf90_inquire_attribute(ncid, NF90_GLOBAL, name)

            link_exists = (ncerr == nf90_noerr)
            ncerr = 0
        end function has_attribute

        function has_dataset(ncid, name) result(link_exists)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            integer                  :: varid
            logical                  :: link_exists

            ncerr = nf90_inq_varid(ncid, name, varid)

            link_exists = (ncerr == nf90_noerr)

            ncerr = nf90_inquire_variable(ncid, varid)

            link_exists = (link_exists .and. (ncerr == nf90_noerr))
            ncerr = 0
        end function has_dataset

        subroutine read_netcdf_dataset_1d_integer(ncid, name, buffer, start, cnt)
            integer,           intent(in)  :: ncid
            character(*),      intent(in)  :: name
            integer,           intent(out) :: buffer(:)
            integer, optional, intent(in)  :: start(:)
            integer, optional, intent(in)  :: cnt(:)
            integer                        :: varid

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading dataset id failed.")

            ncerr = nf90_get_var(ncid=ncid, varid=varid, values=buffer, &
                                 start=start, count=cnt)
        end subroutine read_netcdf_dataset_1d_integer

        subroutine read_netcdf_dataset_1d(ncid, name, buffer, start, cnt)
            integer,           intent(in)  :: ncid
            character(*),      intent(in)  :: name
            double precision,  intent(out) :: buffer(:)
            integer, optional, intent(in)  :: start(:)
            integer, optional, intent(in)  :: cnt(:)
            integer                        :: varid

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading dataset id failed.")

            ncerr = nf90_get_var(ncid=ncid, varid=varid, values=buffer, &
                                 start=start, count=cnt)
        end subroutine read_netcdf_dataset_1d

        subroutine read_netcdf_dataset_2d(ncid, name, buffer, start, cnt)
            integer,           intent(in)  :: ncid
            character(*),      intent(in)  :: name
            double precision,  intent(out) :: buffer(:, :)
            integer, optional, intent(in)  :: start(:)
            integer, optional, intent(in)  :: cnt(:)
            integer                        :: varid, map(2)
            double precision, allocatable  :: values(:, :)

            map = shape(buffer)

            allocate(values(map(2), map(1)))

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading dataset id failed.")

            ncerr = nf90_get_var(ncid=ncid, varid=varid, values=values, &
                                 start=start, count=cnt)

            buffer = transpose(values)

            deallocate(values)
        end subroutine read_netcdf_dataset_2d

        subroutine read_netcdf_dataset_3d(ncid, name, buffer, start, cnt)
            integer,           intent(in)  :: ncid
            character(*),      intent(in)  :: name
            double precision,  intent(out) :: buffer(:, :, :)
            integer, optional, intent(in)  :: start(:)
            integer, optional, intent(in)  :: cnt(:)
            integer                        :: varid, map(3)
            double precision, allocatable  :: values(:, :, :)

            map = shape(buffer)

            allocate(values(map(3), map(2), map(1)))

            ncerr = nf90_inq_varid(ncid, name, varid)
            call check_netcdf_error("Reading dataset id failed.")

            ncerr = nf90_get_var(ncid=ncid, varid=varid, values=values, &
                                 start=start, count=cnt)

            buffer = reshape(values, shape=(/map(1), map(2), map(3)/), order=(/3, 2, 1/))
            deallocate(values)
        end subroutine read_netcdf_dataset_3d

        subroutine read_netcdf_attrib_integer(ncid, name, val)
            integer,       intent(in)     :: ncid
            character(*),  intent(in)     :: name
            integer,       intent(out)    :: val

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, name, val)
            call check_netcdf_error("Reading attribute '" // name // "' failed.")
        end subroutine read_netcdf_attrib_integer

        subroutine read_netcdf_attrib_double(ncid, name, val)
            integer,          intent(in)     :: ncid
            character(*),     intent(in)     :: name
            double precision, intent(out)    :: val

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, name, val)
            call check_netcdf_error("Reading attribute '" // name // "' failed.")

        end subroutine read_netcdf_attrib_double

        subroutine read_netcdf_attrib_character(ncid, name, val)
            integer,      intent(in)     :: ncid
            character(*), intent(in)     :: name
            character(*), intent(out)    :: val

            ncerr = nf90_get_att(ncid, NF90_GLOBAL, name, val)
            call check_netcdf_error("Reading attribute '" // name // "' failed.")

        end subroutine read_netcdf_attrib_character

        subroutine read_netcdf_attrib_logical(ncid, name, val)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            logical,      intent(out) :: val
            character(5)              :: char_val

            call read_netcdf_attrib_character(ncid, name, char_val)

            val = (char_val == 'true')

        end subroutine read_netcdf_attrib_logical

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

        subroutine read_netcdf_attrib_default_double(ncid, name, val)
            integer,          intent(in)   :: ncid
            character(*),     intent(in)    :: name
            double precision, intent(inout) :: val

            if (has_attribute(ncid, name)) then
                call read_netcdf_attrib_double(ncid, name, val)
#ifdef ENABLE_VERBOSE
                print *, "Found float attribute '" // name // "'."
            else
                print *, "WARNING: Using default value of '" // name // "'."
#endif
            endif

        end subroutine read_netcdf_attrib_default_double

        subroutine read_netcdf_attrib_default_integer(ncid, name, val)
            integer,      intent(in)    :: ncid
            character(*), intent(in)    :: name
            integer,      intent(inout) :: val

            if (has_attribute(ncid, name)) then
                call read_netcdf_attrib_integer(ncid, name, val)
#ifdef ENABLE_VERBOSE
                print *, "Found integer attribute '" // name // "'."
            else
                print *, "WARNING: Using default value of '" // name // "'."
#endif
            endif

        end subroutine read_netcdf_attrib_default_integer


        subroutine read_netcdf_attrib_default_logical(ncid, name, val)
            integer,      intent(in)    :: ncid
            character(*), intent(in)    :: name
            logical,      intent(inout) :: val

            if (has_attribute(ncid, name)) then
                call read_netcdf_attrib_logical(ncid, name, val)
#ifdef ENABLE_VERBOSE
                print *, "Found boolean attribute '" // name // "'."
            else
                print *, "WARNING: Using default value of '" // name // "'."
#endif
            endif

        end subroutine read_netcdf_attrib_default_logical

end module netcdf_reader
