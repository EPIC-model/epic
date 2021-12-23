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
    integer            :: ix, iy, ncid, var_id1 = -1, var_id2 = -1, dimids(2)
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

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    call define_netcdf_dimensions(ncid, (/nx, ny/), dimids)

    call define_netcdf_dataset(ncid, 'x_velocity', 'm/s', dimids, var_id1)
    call define_netcdf_dataset(ncid, 'y_velocity', 'm/s', dimids, var_id2)

    call close_definition(ncid)

    call write_netcdf_dataset(ncid, var_id1, dset)

    dset = 10 + dset

    call write_netcdf_dataset(ncid, var_id2, dset)

    call close_netcdf_file(ncid)

    passed = (ncerr == 0)

    call delete_netcdf_file(ncfname='nctest.nc')

    call print_result_logical('Test netCDF write 2D datasets', passed)

end program test_netcdf_dataset_2d
