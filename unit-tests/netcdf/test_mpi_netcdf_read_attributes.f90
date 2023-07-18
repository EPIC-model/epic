! =============================================================================
!                       Test MPI netCDF read attributes
!
!               This unit test checks reading attributes like file type
!               and time to perform restarts in parallel.
! =============================================================================
program test_mpi_netcdf_read_attributes
    use unit_test
    use netcdf_writer
    use netcdf_reader
    use mpi_communicator
    implicit none

    integer, parameter :: nx = 5, ny = 10, nz = 20
    integer            :: ncid, t_dim_id, t_axis_id, n
    logical            :: passed = .true.
    character(len=16)  :: file_type
    double precision   :: t

    call mpi_comm_initialise

    !
    !
    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call write_netcdf_attribute(ncid=ncid, name='file_type', val='fields')
    call write_netcdf_box(ncid, (/0.0d0, 0.0d0, 0.0d0/), &
                                (/1.0d0, 1.0d0, 1.0d0/), &
                                (/nx, ny, nz/))

    call define_netcdf_dimension(ncid=ncid,                         &
                                 name='t',                          &
                                 dimsize=NF90_UNLIMITED,            &
                                 dimid=t_dim_id)
    call define_netcdf_dataset(                                     &
        ncid=ncid,                                                  &
        name='t',                                                   &
        long_name='time ',                                          &
        std_name='time',                                            &
        unit='s',                                                   &
        dtype=NF90_DOUBLE,                                          &
        dimids=(/t_dim_id/),                                        &
        varid=t_axis_id)
    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)
    do n = 1, 10
        call write_netcdf_scalar(ncid, t_axis_id, dble(n), n)
    enddo

    call close_netcdf_file(ncid)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_NOWRITE, &
                          ncid=ncid)

    call get_file_type(ncid, file_type)

    call get_time(ncid, t)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_NOWRITE, &
                          ncid=ncid)

    call get_file_type(ncid, file_type)

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0) .and. (file_type == 'fields') .and. (t == 10.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    !
    !
    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call write_netcdf_attribute(ncid=ncid, name='file_type', val='field_stats')
    call write_netcdf_box(ncid, (/0.0d0, 0.0d0, 0.0d0/), &
                                (/1.0d0, 1.0d0, 1.0d0/), &
                                (/nx, ny, nz/))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_NOWRITE, &
                          ncid=ncid)

    call get_file_type(ncid, file_type)

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0) .and. (file_type == 'field_stats'))

    call delete_netcdf_file(ncfname='nctest.nc')

    !
    !
    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    call write_netcdf_attribute(ncid=ncid, name='file_type', val='parcels')
    call write_netcdf_box(ncid, (/0.0d0, 0.0d0, 0.0d0/), &
                                (/1.0d0, 1.0d0, 1.0d0/), &
                                (/nx, ny, nz/))

    call write_netcdf_scalar(ncid, t_axis_id, 12.0, 1)

    call close_netcdf_file(ncid)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_NOWRITE, &
                          ncid=ncid)

    call get_file_type(ncid, file_type)

    call get_time(ncid, t)

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0) .and. (file_type == 'parcels') .and. (t == 12.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    !
    !
    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call write_netcdf_attribute(ncid=ncid, name='file_type', val='parcel_stats')
    call write_netcdf_box(ncid, (/0.0d0, 0.0d0, 0.0d0/), &
                                (/1.0d0, 1.0d0, 1.0d0/), &
                                (/nx, ny, nz/))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_NOWRITE, &
                          ncid=ncid)

    call get_file_type(ncid, file_type)

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0) .and. (file_type == 'parcel_stats'))

    call delete_netcdf_file(ncfname='nctest.nc')

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI netCDF read attributes', passed)
    endif

    call mpi_comm_finalise

end program test_mpi_netcdf_read_attributes
