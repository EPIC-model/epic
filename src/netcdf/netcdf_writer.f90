module netcdf_writer
    use netcdf_utils
    implicit none

    private :: write_netcdf_dataset_1d, &
               write_netcdf_dataset_2d, &
               write_netcdf_dataset_3d

    interface write_netcdf_dataset
        module procedure write_netcdf_dataset_1d
        module procedure write_netcdf_dataset_2d
        module procedure write_netcdf_dataset_3d
    end interface write_netcdf_dataset

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

        subroutine write_netcdf_scalar(ncid, varid, data)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data

            ! write data
            ncerr = nf90_put_var(ncid, varid, data)

            call check_netcdf_error("Failed to write scalar.")

        end subroutine write_netcdf_scalar

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
            ncerr = nf90_put_var(ncid, varid, data)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_3d

end module netcdf_writer
