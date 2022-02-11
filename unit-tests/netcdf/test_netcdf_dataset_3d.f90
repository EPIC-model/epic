! =============================================================================
!                       Test netCDF 3D dataset
!
!               This unit test checks to write 3D datasets to netCDF.
! =============================================================================
program test_netcdf_dataset_3d
    use unit_test
    use netcdf_writer
    implicit none

    integer, parameter :: nx = 5, ny = 10, nz = 2
    integer            :: ix, iy, iz, ncid, dimids(3)
    integer            :: var_id1 = -1, var_id2 = -1, var_id3
    double precision   :: dset(nz, ny, nx)
    logical            :: passed

    do ix = 1, nx
        do iy = 1, ny
            do iz = 1, nz
                dset(iz, iy, ix) = iz + (iy-1) * nz + (ix-1) * ny * nz
            enddo
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    call define_netcdf_dimension(ncid, "x", nx, dimids(3))
    call define_netcdf_dimension(ncid, "y", ny, dimids(2))
    call define_netcdf_dimension(ncid, "z", nz, dimids(1))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id1)
    call define_netcdf_dataset(ncid, 'y_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id2)
    call define_netcdf_dataset(ncid, 'nparg', '', '', '-', NF90_INT, dimids, var_id3)

    call close_definition(ncid)

    call write_netcdf_dataset(ncid, var_id1, dset)

    dset = 1.5d0 + dset

    call write_netcdf_dataset(ncid, var_id2, dset)

    call write_netcdf_dataset(ncid, var_id3, dset)

    call close_netcdf_file(ncid)

    passed = (ncerr == 0)

!     call delete_netcdf_file(ncfname='nctest.nc')

    call print_result_logical('Test netCDF write 3D datasets', passed)

end program test_netcdf_dataset_3d
