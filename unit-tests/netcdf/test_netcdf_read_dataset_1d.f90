! =============================================================================
!                       Test netCDF reading 1D datasets
!
!               This unit test checks to read 1D datasets. It checks
!               reading in a dataset into a bigger array.
! =============================================================================
program test_netcdf_read_dataset_1d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    implicit none

    integer, parameter :: n_parcels = 100
    integer            :: n, ncid, dimid
    integer            :: var_id = -1, cnt(1), start(1)
    double precision   :: wdset(n_parcels), rdset(n_parcels + 100), error
    logical            :: passed

    do n = 1, n_parcels
        wdset(n) = dble(n)
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    call define_netcdf_dimension(ncid, "n_parcels", n_parcels, dimid)

    call define_netcdf_dataset(ncid, 'wdset', '', '', '', NF90_DOUBLE, (/dimid/), var_id)

    call close_definition(ncid)

    ! write data
    cnt   = (/ n_parcels /)
    start = (/ 1 /)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)
    call write_netcdf_dataset(ncid, var_id, wdset, start, cnt)


    call close_netcdf_file(ncid)

    ! read data
    call open_netcdf_file(ncfname='nctest.nc',        &
                            access_flag=NF90_NOWRITE, &
                            ncid=ncid)

    call read_netcdf_dataset(ncid, 'wdset', rdset(1:n_parcels), start, cnt)

    call close_netcdf_file(ncid)

    error = sum(wdset - rdset(1:n_parcels))

    passed = ((ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    call print_result_logical('Test netCDF read 1D dataset', passed)

end program test_netcdf_read_dataset_1d
