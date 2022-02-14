module netcdf_writer
    use netcdf_utils
    implicit none

    private :: write_netcdf_dataset_1d,             &
               write_netcdf_dataset_2d,             &
               write_netcdf_dataset_3d,             &
               write_netcdf_scalar_integer,         &
               write_netcdf_scalar_double,          &
               write_netcdf_attribute_integer,      &
               write_netcdf_attribute_double,       &
               write_netcdf_attribute_character,    &
               write_netcdf_attribute_logical

    interface write_netcdf_scalar
        module procedure write_netcdf_scalar_integer
        module procedure write_netcdf_scalar_double
    end interface write_netcdf_scalar

    interface write_netcdf_dataset
        module procedure write_netcdf_dataset_1d
        module procedure write_netcdf_dataset_2d
        module procedure write_netcdf_dataset_3d
    end interface write_netcdf_dataset

    interface write_netcdf_attribute
        module procedure write_netcdf_attribute_integer
        module procedure write_netcdf_attribute_double
        module procedure write_netcdf_attribute_character
        module procedure write_netcdf_attribute_logical
    end interface write_netcdf_attribute

    contains

        subroutine define_netcdf_dimension(ncid, name, dimsize, dimid)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            integer,      intent(in)  :: dimsize
            integer,      intent(out) :: dimid

            ncerr = nf90_def_dim(ncid, name, dimsize, dimid)
            call check_netcdf_error("Failed to define" // name // "dimension.")
        end subroutine define_netcdf_dimension

        subroutine define_netcdf_dataset(ncid, name, long_name, std_name, unit, dtype, dimids, varid)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            character(*), intent(in)  :: long_name
            character(*), intent(in)  :: std_name
            character(*), intent(in)  :: unit
            integer,      intent(in)  :: dtype ! NF90_DOUBLE or NF90_INT
            integer,      intent(in)  :: dimids(:)
            integer,      intent(out) :: varid

            ! define the variable
            ncerr = nf90_def_var(ncid, name, dtype, dimids, varid)

            call check_netcdf_error("Failed to define the dataset.")

            ncerr = nf90_put_att(ncid, varid, "units", unit)

            call check_netcdf_error("Failed to define the dataset unit.")

            if (len(long_name) > 0) then
                ncerr = nf90_put_att(ncid, varid, "long_name", long_name)
                call check_netcdf_error("Failed to define the dataset long name.")
            endif

            if (len(std_name) > 0) then
                ncerr = nf90_put_att(ncid, varid, "standard_name", std_name)
                call check_netcdf_error("Failed to define the dataset standard name.")
            endif

        end subroutine define_netcdf_dataset

        subroutine close_definition(ncid)
            integer, intent(in) :: ncid

            ! tell netCDF metadata is defined
            ncerr = nf90_enddef(ncid)

            call check_netcdf_error("Failed to define metadata.")
        end subroutine close_definition

        subroutine write_netcdf_attribute_integer(ncid, name, val)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            integer,      intent(in) :: val

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name=name, values=val)
            call check_netcdf_error("Failed to define '" // name // "' attribute.")
        end subroutine write_netcdf_attribute_integer

        subroutine write_netcdf_attribute_double(ncid, name, val)
            integer,          intent(in) :: ncid
            character(*),     intent(in) :: name
            double precision, intent(in) :: val

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name=name, values=val)
            call check_netcdf_error("Failed to define '" // name // "' attribute.")
        end subroutine write_netcdf_attribute_double

        subroutine write_netcdf_attribute_character(ncid, name, val)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            character(*), intent(in) :: val

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name=name, values=val)
            call check_netcdf_error("Failed to define '" // name // "' attribute.")
        end subroutine write_netcdf_attribute_character

        subroutine write_netcdf_attribute_logical(ncid, name, val)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: name
            logical,      intent(in) :: val
            integer                  :: int_val

            int_val = 0
            if (val) then
                int_val = 1
            endif

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name=name, values=int_val)
            call check_netcdf_error("Failed to define '" // name // "' attribute.")

        end subroutine write_netcdf_attribute_logical

        subroutine write_netcdf_scalar_integer(ncid, varid, data, start)
            integer, intent(in) :: ncid
            integer, intent(in) :: varid
            integer, intent(in) :: data
            integer, intent(in) :: start

            ! write data
            ncerr = nf90_put_var(ncid, varid, (/data/), &
                                 start=(/start/), count=(/1/))

            call check_netcdf_error("Failed to write scalar.")

        end subroutine write_netcdf_scalar_integer

        subroutine write_netcdf_scalar_double(ncid, varid, data, start)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data
            integer,           intent(in) :: start

            ! write data
            ncerr = nf90_put_var(ncid, varid, (/data/), &
                                 start=(/start/), count=(/1/))

            call check_netcdf_error("Failed to write scalar.")

        end subroutine write_netcdf_scalar_double

        subroutine write_netcdf_dataset_1d(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            ! write data
            ncerr = nf90_put_var(ncid, varid, data, &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_1d

        subroutine write_netcdf_dataset_2d(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            ! write data
            ncerr = nf90_put_var(ncid, varid, data, &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_2d

        subroutine write_netcdf_dataset_3d(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:, :, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            ! write data
            ncerr = nf90_put_var(ncid, varid, data, &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_3d

        subroutine write_netcdf_timestamp(ncid)
            integer,          intent(in) :: ncid
            character(len=8)             :: cdate, tmp2
            character(len=10)            :: ctime, tmp1
            character(len=5)             :: czone
            character(len=9)             :: tmp3

            ! 15 June 2021
            ! https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005fTIME.html
            call date_and_time(date=cdate, time=ctime, zone=czone)
            ! 15 June 2021
            ! https://stackoverflow.com/questions/13755762/access-character-at-specific-index-in-a-string-in-fortran
            tmp1 = cdate(1:4) // '/' // cdate(5:6) // '/' // cdate(7:8)
            tmp2 = ctime(1:2) // ':' // ctime(3:4) // ':' // ctime(5:6)
            tmp3 = 'UTC' // czone(1:3) // ':' // czone(4:5)
            call write_netcdf_attribute_character(ncid, "creation_date", tmp1)
            call write_netcdf_attribute_character(ncid, "creation_time", tmp2)
            call write_netcdf_attribute_character(ncid, "creation_zone", tmp3)
        end subroutine write_netcdf_timestamp

        subroutine write_netcdf_box(ncid, origin, extent, ncells)
            integer,          intent(in) :: ncid
            double precision, intent(in) :: origin(:), extent(:)
            integer,          intent(in) :: ncells(:)

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name="ncells", values=ncells)
            call check_netcdf_error("Failed to define 'ncells' global attribute.")

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name="extent", values=extent)
            call check_netcdf_error("Failed to define 'extent' global attribute.")

            ncerr = nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, name="origin", values=origin)
            call check_netcdf_error("Failed to define 'origin' global attribute.")
        end subroutine write_netcdf_box

end module netcdf_writer
