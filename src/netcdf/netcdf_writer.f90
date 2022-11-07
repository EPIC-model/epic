! Note: In parallel we need to set NF90_COLLECTIVE for variables if one dimension is NF90_UNLIMITED
module netcdf_writer
    use netcdf_utils
    implicit none

    private :: write_netcdf_dataset_1d_double,      &
               write_netcdf_dataset_2d_double,      &
               write_netcdf_dataset_3d_double,      &
               write_netcdf_dataset_1d_integer,     &
               write_netcdf_dataset_2d_integer,     &
               write_netcdf_dataset_3d_integer,     &
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
        module procedure write_netcdf_dataset_1d_double
        module procedure write_netcdf_dataset_2d_double
        module procedure write_netcdf_dataset_3d_double
        module procedure write_netcdf_dataset_1d_integer
        module procedure write_netcdf_dataset_2d_integer
        module procedure write_netcdf_dataset_3d_integer
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
            character(5)             :: char_val

            char_val = 'false'
            if (val) then
                char_val = 'true'
            endif

            call write_netcdf_attribute_character(ncid, name, trim(char_val))

        end subroutine write_netcdf_attribute_logical

        subroutine write_netcdf_scalar_integer(ncid, varid, data, start)
            integer, intent(in) :: ncid
            integer, intent(in) :: varid
            integer, intent(in) :: data
            integer, intent(in) :: start

            if (mpi_size > 1) then
                ! 18 April 2022
                ! see also https://github.com/Unidata/netcdf-fortran/blob/main/examples/F90/simple_xy_par_wr2.f90
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

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

            if (mpi_size > 1) then
                ! 18 April 2022
                ! see also https://github.com/Unidata/netcdf-fortran/blob/main/examples/F90/simple_xy_par_wr2.f90
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! write data
            ncerr = nf90_put_var(ncid, varid, (/data/), &
                                 start=(/start/), count=(/1/))

            call check_netcdf_error("Failed to write scalar.")

        end subroutine write_netcdf_scalar_double

        subroutine write_netcdf_dataset_1d_double(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            if (mpi_size > 1) then
                ! 18 April 2022
                ! see also https://github.com/Unidata/netcdf-fortran/blob/main/examples/F90/simple_xy_par_wr2.f90
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! write data
            ncerr = nf90_put_var(ncid, varid, data, &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_1d_double

        subroutine write_netcdf_dataset_2d_double(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            ! transpose(data) == reshape(data, shape=(/map(2), map(1)/), order=(/2, 1/)
            ! with map = shape(data)

            if (mpi_size > 1) then
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! write data
            ncerr = nf90_put_var(ncid, varid, transpose(data),  &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_2d_double

        subroutine write_netcdf_dataset_3d_double(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            double precision,  intent(in) :: data(:, :, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)
            integer                       :: map(3)

            if (mpi_size > 1) then
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            map = shape(data)

            ! write data
            ncerr = nf90_put_var(ncid, varid,                                    &
                                 reshape(data, shape=(/map(3), map(2), map(1)/), &
                                         order=(/3, 2, 1/)),                     &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_3d_double

        subroutine write_netcdf_dataset_1d_integer(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            integer,           intent(in) :: data(:)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            if (mpi_size > 1) then
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! write data
            ncerr = nf90_put_var(ncid, varid, data, &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_1d_integer

        subroutine write_netcdf_dataset_2d_integer(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            integer,           intent(in) :: data(:, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)

            if (mpi_size > 1) then
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! transpose(data) == reshape(data, shape=(/map(2), map(1)/), order=(/2, 1/)
            ! with map = shape(data)

            ! write data
            ncerr = nf90_put_var(ncid, varid, transpose(data),  &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_2d_integer

        subroutine write_netcdf_dataset_3d_integer(ncid, varid, data, start, cnt)
            integer,           intent(in) :: ncid
            integer,           intent(in) :: varid
            integer,           intent(in) :: data(:, :, :)
            integer, optional, intent(in) :: start(:)
            integer, optional, intent(in) :: cnt(:)
            integer                       :: map(3)

            map = shape(data)

            if (mpi_size > 1) then
                ncerr = nf90_var_par_access(ncid, varid, NF90_COLLECTIVE)
                call check_netcdf_error("Failed to set collective.")
            endif

            ! write data
            ncerr = nf90_put_var(ncid, varid,                                    &
                                 reshape(data, shape=(/map(3), map(2), map(1)/), &
                                         order=(/3, 2, 1/)),                     &
                                 start=start, count = cnt)

            call check_netcdf_error("Failed to write dataset.")

        end subroutine write_netcdf_dataset_3d_integer

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

        subroutine write_netcdf_axis_2d(ncid, dimids, origin, dx, ngps)
            integer,          intent(in) :: ncid
            double precision, intent(in) :: origin(2), dx(2)
            integer,          intent(in) :: dimids(2), ngps(2)
            integer                      :: i
            double precision             :: x_axis(0:ngps(1)-1), z_axis(0:ngps(2)-1)

            do i = 0, ngps(1)-1
                x_axis(i) = origin(1) + dble(i) * dx(1)
            enddo

            do i = 0, ngps(2)-1
                z_axis(i) = origin(2) + dble(i) * dx(2)
            enddo

            call write_netcdf_dataset(ncid, dimids(1), x_axis)
            call write_netcdf_dataset(ncid, dimids(2), z_axis)
        end subroutine write_netcdf_axis_2d

        subroutine write_netcdf_axis_3d(ncid, dimids, origin, dx, ngps)
            integer,          intent(in) :: ncid
            double precision, intent(in) :: origin(3), dx(3)
            integer,          intent(in) :: dimids(3), ngps(3)
            integer                      :: i
            double precision             :: x_axis(0:ngps(1)-1)
            double precision             :: y_axis(0:ngps(2)-1)
            double precision             :: z_axis(0:ngps(3)-1)


            do i = 0, ngps(1)-1
                x_axis(i) = origin(1) + dble(i) * dx(1)
            enddo

            do i = 0, ngps(2)-1
                y_axis(i) = origin(2) + dble(i) * dx(2)
            enddo

            do i = 0, ngps(3)-1
                z_axis(i) = origin(3) + dble(i) * dx(3)
            enddo

            call write_netcdf_dataset(ncid, dimids(1), x_axis)
            call write_netcdf_dataset(ncid, dimids(2), y_axis)
            call write_netcdf_dataset(ncid, dimids(3), z_axis)
        end subroutine write_netcdf_axis_3d

        subroutine write_netcdf_info(ncid, version_tag, file_type, cf_version)
            integer,      intent(in) :: ncid
            character(*), intent(in) :: version_tag
            character(*), intent(in) :: file_type
            character(*), intent(in) :: cf_version
            call write_netcdf_attribute(ncid=ncid, name=version_name, val=version_tag)
            call write_netcdf_attribute(ncid=ncid, name='file_type', val=file_type)
            call write_netcdf_attribute(ncid=ncid, name='Conventions', val=cf_version)
            call write_netcdf_timestamp(ncid)
        end subroutine write_netcdf_info

        subroutine define_netcdf_spatial_dimensions_2d(ncid, ngps, dimids, axids)
            integer, intent(in)  :: ncid
            integer, intent(in)  :: ngps(2)
            integer, intent(out) :: dimids(2), axids(2)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                 &
                                         name=netcdf_dims(1),       &
                                         dimsize=ngps(1),           &
                                         dimid=dimids(1))

            call define_netcdf_dimension(ncid=ncid,                 &
                                         name=netcdf_dims(2),       &
                                         dimsize=ngps(2),           &
                                         dimid=dimids(2))

            call define_netcdf_dataset(                                                 &
                ncid=ncid,                                                              &
                name=netcdf_dims(1),                                                    &
                long_name=netcdf_dims(1)//'-coordinate in projected coordinate system', &
                std_name='projection_'//netcdf_dims(1)//'_coordinate',                  &
                unit='m',                                                               &
                dtype=NF90_DOUBLE,                                                      &
                dimids=(/dimids(1)/),                                                   &
                varid=axids(1))

            ncerr = nf90_put_att(ncid, axids(1), "axis", netcdf_axes(1))
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                         &
                ncid=ncid,                                                      &
                name=netcdf_dims(2),                                            &
                long_name='height coordinate in projected coordinate system',   &
                std_name='height',                                              &
                unit='m',                                                       &
                dtype=NF90_DOUBLE,                                              &
                dimids=(/dimids(2)/),                                           &
                varid=axids(2))

            ncerr = nf90_put_att(ncid, axids(2), "axis", netcdf_axes(2))
            call check_netcdf_error("Failed to add axis attribute.")

        end subroutine define_netcdf_spatial_dimensions_2d

        subroutine define_netcdf_spatial_dimensions_3d(ncid, ngps, dimids, axids)
            integer, intent(in)  :: ncid
            integer, intent(in)  :: ngps(3)
            integer, intent(out) :: dimids(3), axids(3)

            call define_netcdf_dimension(ncid=ncid,           &
                                         name=netcdf_dims(1), &
                                         dimsize=ngps(1),     &
                                         dimid=dimids(1))

            call define_netcdf_dimension(ncid=ncid,           &
                                         name=netcdf_dims(2), &
                                         dimsize=ngps(2),     &
                                         dimid=dimids(2))

            call define_netcdf_dimension(ncid=ncid,           &
                                         name=netcdf_dims(3), &
                                         dimsize=ngps(3),     &
                                         dimid=dimids(3))

            call define_netcdf_dataset(                                                 &
                ncid=ncid,                                                              &
                name=netcdf_dims(1),                                                    &
                long_name=netcdf_dims(1)//'-coordinate in projected coordinate system', &
                std_name='projection_'//netcdf_dims(1)//'_coordinate',                  &
                unit='m',                                                               &
                dtype=NF90_DOUBLE,                                                      &
                dimids=(/dimids(1)/),                                                   &
                varid=axids(1))

            ncerr = nf90_put_att(ncid, axids(1), "axis", netcdf_axes(1))
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                                 &
                ncid=ncid,                                                              &
                name=netcdf_dims(2),                                                    &
                long_name=netcdf_dims(2)//'-coordinate in projected coordinate system', &
                std_name='projection_'//netcdf_dims(2)//'_coordinate',                  &
                unit='m',                                                               &
                dtype=NF90_DOUBLE,                                                      &
                dimids=(/dimids(2)/),                                                   &
                varid=axids(2))

            ncerr = nf90_put_att(ncid, axids(2), "axis", netcdf_axes(2))
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                             &
                ncid=ncid,                                                          &
                name=netcdf_dims(3),                                                &
                long_name='height coordinate in projected coordinate system',       &
                std_name='height',                                                  &
                unit='m',                                                           &
                dtype=NF90_DOUBLE,                                                  &
                dimids=(/dimids(3)/),                                               &
                varid=axids(3))

            ncerr = nf90_put_att(ncid, axids(3), "axis", netcdf_axes(3))
            call check_netcdf_error("Failed to add axis attribute.")

        end subroutine define_netcdf_spatial_dimensions_3d

        subroutine define_netcdf_temporal_dimension(ncid, dimid, axid)
            integer, intent(in)  :: ncid
            integer, intent(out) :: dimid, axid

            call define_netcdf_dimension(ncid=ncid,                 &
                                         name=netcdf_dims(4),       &
                                         dimsize=NF90_UNLIMITED,    &
                                         dimid=dimid)

            call define_netcdf_dataset(             &
                ncid=ncid,                          &
                name=netcdf_dims(4),                &
                long_name='time ',                  &
                std_name='time',                    &
                unit='seconds since 1970-01-01',    &
                dtype=NF90_DOUBLE,                  &
                dimids=(/dimid/),                   &
                varid=axid)

            ncerr = nf90_put_att(ncid, axid, "axis", netcdf_axes(4))
            call check_netcdf_error("Failed to add axis attribute.")

            ncerr = nf90_put_att(ncid, axid, "calendar", &
                                 'proleptic_gregorian')
            call check_netcdf_error("Failed to add calendear attribute.")
        end subroutine define_netcdf_temporal_dimension

end module netcdf_writer
