! =============================================================================
!                       Test MPI netCDF 1D dataset
!
!       This unit test checks to write 2D datasets in parallel to netCDF.
! =============================================================================
program test_mpi_netcdf_dataset_2d
    use unit_test
    use netcdf_writer
    use mpi_communicator
    implicit none

    integer, parameter :: n_local_parcels = 5
    integer            :: n_total_parcels
    integer            :: n, ncid, dimids
    integer            :: var_id = -1
    integer            :: dset(n_local_parcels)
    logical            :: passed = .true.
    integer            :: start(1), cnt(1)

    call mpi_comm_initialise

    do n = 1, n_local_parcels
        dset(n) = n + n_local_parcels * comm%rank
    enddo

    n_total_parcels = comm%size * n_local_parcels

    call create_netcdf_file(ncfname='nctest.nc',                    &
                            overwrite=.true.,                       &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "n_parcels", n_total_parcels, dimids)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'dummy', '', '', '1', NF90_INT, (/dimids/), var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',          &
                          access_flag=NF90_WRITE,       &
                          ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    start(1) = comm%rank * n_local_parcels + 1
    cnt(1) = n_local_parcels

    call write_netcdf_dataset(ncid, var_id, dset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI netCDF write 1D dataset', passed)
    endif

    call mpi_comm_finalise

end program test_mpi_netcdf_dataset_2d
