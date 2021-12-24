! =============================================================================
!                       Test netCDF time series
!
!               This unit test checks to write multiple time steps.
! =============================================================================
program test_netcdf_time_series_2
    use unit_test
    use netcdf_writer
    implicit none

    integer, parameter :: nx = 5, ny = 10, nt = 3
    integer            :: ix, iy, ncid, dimids(3), t
    integer            :: var_id = -1, cnt(3), start(3)
    double precision   :: dset(ny, nx)
    logical            :: passed

    do ix = 1, nx
        do iy = 1, ny
            dset(iy, ix) = iy + (ix-1) * ny
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call define_netcdf_dimension(ncid, "y", ny, dimids(1))
    call define_netcdf_dimension(ncid, "x", nx, dimids(2))
    call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, dimids(3))

    call define_netcdf_dataset(ncid, 'x_velocity', 'm/s', NF90_DOUBLE, dimids, var_id)

    call close_definition(ncid)


    do t = 1, nt

        cnt   = (/ ny, nx, 1 /)
        start = (/ 1,  1,  t /)

        call open_netcdf_file(ncfname='nctest.nc',    &
                            access_flag=NF90_WRITE, &
                            ncid=ncid)
        call write_netcdf_dataset(ncid, var_id, dset, start, cnt)


        call close_netcdf_file(ncid)

        dset = 1.0d0 + dset
    enddo

    passed = (ncerr == 0)

    call delete_netcdf_file(ncfname='nctest.nc')

    call print_result_logical('Test netCDF write time series 2', passed)

end program test_netcdf_time_series_2
