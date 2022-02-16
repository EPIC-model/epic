! =============================================================================
!                       Test netCDF reading 3D datasets
!
!               This unit test checks to read 3D datasets.
! =============================================================================
program test_netcdf_read_dataset_3d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    implicit none

    integer, parameter :: nx = 5, ny = 10, nz = 15, nt = 3
    integer            :: ix, iy, iz, ncid, dimids(4), t
    integer            :: var_id = -1, cnt(4), start(4)
    double precision   :: wdset(nz, ny, nx), rdset(nz, ny, nx), error
    logical            :: passed = .true.

    do ix = 1, nx
        do iy = 1, ny
            do iz = 1, nz
                wdset(iz, iy, ix) = iy + (ix-1) * ny + 3 * iz
            enddo
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

    call define_netcdf_dimension(ncid, "z", nz, dimids(3))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, dimids(4))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    ! write data

    do t = 1, nt

        cnt   = (/ nx, ny, nz, 1 /)
        start = (/ 1,  1,  1,  t /)

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

    cnt   = (/nx, ny, nz, 1/)
    start = (/1,  1,  1,  nt/)
    call read_netcdf_dataset(ncid, 'x_velocity', rdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    error = sum(wdset - rdset)

    passed = (passed .and. (ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    call print_result_logical('Test netCDF read 3D dataset', passed)

end program test_netcdf_read_dataset_3d
