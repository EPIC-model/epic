! =============================================================================
!                       Test netCDF reading 2D datasets
!
!               This unit test checks to read 2D datasets.
! =============================================================================
program test_netcdf_read_dataset_2d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    implicit none

    integer, parameter :: nx = 5, ny = 10, nt = 3
    integer            :: ix, iy, ncid, dimids(3), t
    integer            :: var_id = -1, cnt(3), start(3)
    double precision   :: wdset(ny, nx), rdset(ny, nx), error
    logical            :: passed = .true.

    do ix = 1, nx
        do iy = 1, ny
            wdset(iy, ix) = iy + (ix-1) * ny
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

    call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, dimids(3))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    ! write data

    do t = 1, nt

        cnt   = (/ nx, ny, 1 /)
        start = (/ 1,  1,  t /)

        call open_netcdf_file(ncfname='nctest.nc',    &
                            access_flag=NF90_WRITE, &
                            ncid=ncid)

        passed = (passed .and. (ncerr == 0))

        call write_netcdf_dataset(ncid, var_id, wdset, start, cnt)

        passed = (passed .and. (ncerr == 0))

        call close_netcdf_file(ncid)

        passed = (passed .and. (ncerr == 0))

        wdset = 1.0d0 + wdset
    enddo

    ! we need to subtract the last increase again
    ! since it was not written
    wdset = wdset - 1.0d0

    ! read last data
    call open_netcdf_file(ncfname='nctest.nc',        &
                            access_flag=NF90_NOWRITE, &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    cnt   = (/nx, ny, 1/)
    start = (/1,  1,  nt/)
    call read_netcdf_dataset(ncid, 'x_velocity', rdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    error = sum(abs(wdset - rdset))

    passed = (passed .and. (ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    call print_result_logical('Test netCDF read 2D dataset', passed)

end program test_netcdf_read_dataset_2d
