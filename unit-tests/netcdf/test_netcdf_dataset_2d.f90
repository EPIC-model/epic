! =============================================================================
!                       Test netCDF 2D dataset
!
!               This unit test checks to write 2D datasets to netCDF.
! =============================================================================
program test_netcdf_dataset_2d
    use unit_test
    use netcdf_writer
    implicit none

    integer, parameter :: nx = 5, ny = 10
    integer            :: ix, iy, ncid, dimids(2), map(2)
    integer            :: var_id1 = -1, var_id2 = -1, var_id3 = -1
    double precision   :: dset(ny, nx)
    logical            :: passed = .true.

    do ix = 1, nx
        do iy = 1, ny
            dset(iy, ix) = iy + (ix-1) * ny
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "x", nx, dimids(1))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "y", ny, dimids(2))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id1)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'y_velocity', '', '',  'm/s', NF90_DOUBLE, dimids, var_id2)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'nparg', '', '',  '-', NF90_INT, dimids, var_id3)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id1, dset)

    passed = (passed .and. (ncerr == 0))

    dset = 1.5d0 + dset

    call write_netcdf_dataset(ncid, var_id2, dset)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id3, dset)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    call print_result_logical('Test netCDF write 2D datasets', passed)

end program test_netcdf_dataset_2d
