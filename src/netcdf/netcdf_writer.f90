module netcdf_writer
    use netcdf_utils
    implicit none

    private :: write_netcdf_dataset_2d, &
               write_netcdf_dataset_3d

    interface write_netcdf_dataset
        module procedure write_netcdf_dataset_2d
        module procedure write_netcdf_dataset_3d
    end interface write_netcdf_dataset

    contains

        subroutine define_netcdf_dimensions(ncid, dims, dimids)
            integer, intent(in)  :: ncid
            integer, intent(in)  :: dims(:)
            integer, intent(out) :: dimids(:)

            ! the fields have ordering: z, y, x; hence we need to
            ! change dimids access

            if (size(dims) == 2) then
                ncerr = nf90_def_dim(ncid, "x", dims(1), dimids(2))
                call check_netcdf_error("Failed to define x dimension.")

                ncerr = nf90_def_dim(ncid, "z", dims(2), dimids(1))
                call check_netcdf_error("Failed to define z dimension.")

            else if (size(dims) == 3) then
                ncerr = nf90_def_dim(ncid, "x", dims(1), dimids(3))
                call check_netcdf_error("Failed to define x dimension.")

                ncerr = nf90_def_dim(ncid, "y", dims(2), dimids(2))
                call check_netcdf_error("Failed to define y dimension.")

                ncerr = nf90_def_dim(ncid, "z", dims(3), dimids(1))
                call check_netcdf_error("Failed to define z dimension.")
            else
                print *, 'Can only define 2D and 3D.'
                stop
            endif
        end subroutine define_netcdf_dimensions

        subroutine define_netcdf_dataset(ncid, name, unit, dtype, dimids, varid)
            integer,      intent(in)  :: ncid
            character(*), intent(in)  :: name
            character(*), intent(in)  :: unit
            integer,      intent(in)  :: dtype ! NF90_DOUBLE or NF90_INT
            integer,      intent(in)  :: dimids(:)
            integer,      intent(out) :: varid

            ! define the variable
            ncerr = nf90_def_var(ncid, name, dtype, dimids, varid)

            call check_netcdf_error("Failed to define the dataset.")

            ncerr = nf90_put_att(ncid, varid, "units", unit)

            call check_netcdf_error("Failed to define the dataset unit.")

        end subroutine define_netcdf_dataset

        subroutine close_definition(ncid)
            integer, intent(in) :: ncid

            ! tell netCDF metadata is defined
            ncerr = nf90_enddef(ncid)

            call check_netcdf_error("Failed to define metadata.")
        end subroutine close_definition

        subroutine write_netcdf_dataset_2d(ncid, varid, data)
            integer, intent(in)          :: ncid
            integer, intent(in)          :: varid
            double precision, intent(in) :: data(:, :)

            ! write data
            ncerr = nf90_put_var(ncid, varid, data)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_2d

        subroutine write_netcdf_dataset_3d(ncid, varid, data)
            integer, intent(in)          :: ncid
            integer, intent(in)          :: varid
            double precision, intent(in) :: data(:, :, :)

            ! write data
            ncerr = nf90_put_var(ncid, varid, data)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_3d

end module netcdf_writer
