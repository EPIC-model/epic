! =============================================================================
!                       Test MPI netCDF reading 1D datasets
!
!               This unit test checks to read 1D datasets in parallel.
! =============================================================================
program test_mpi_netcdf_read_dataset_1d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    use mpi_communicator
    implicit none

    integer, parameter            :: n_local_parcels = 10
    integer                       :: n_parcels
    integer                       :: n, ncid, dimid
    integer                       :: var_id = -1, cnt(1), start(1)
    double precision, allocatable :: wdset(:), rdset(:), error
    logical                       :: passed = .true.
    integer, allocatable          :: recvcounts(:), n_size(:)
    integer                       :: recvbuf

    call mpi_comm_initialise

    n_parcels = n_local_parcels * mpi_size

    allocate(wdset(n_local_parcels))
    allocate(rdset(n_local_parcels))
    allocate(recvcounts(mpi_size))
    allocate(n_size(mpi_size))

    do n = 1, n_local_parcels
        wdset(n) = dble(n_local_parcels * mpi_rank + n)
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "n_parcels", n_parcels, dimid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'wdset', '', '', '', NF90_DOUBLE, (/dimid/), var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))


    ! get start index of writing
    recvcounts = 1
    recvbuf = 0
    n_size = 0
    n_size(mpi_rank+1:mpi_size) = n_local_parcels
    n_size(mpi_rank+1) = 0

    call MPI_Reduce_scatter(n_size, recvbuf, recvcounts, MPI_INT, MPI_SUM, comm_world, mpi_err)

    ! time step to write [step(2) is the time]
    cnt   = (/ n_local_parcels /)
    start = (/ 1 + recvbuf     /)

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id, wdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    ! read data
    call open_netcdf_file(ncfname='nctest.nc',        &
                            access_flag=NF90_NOWRITE, &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call read_netcdf_dataset(ncid, 'wdset', rdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    error = sum(abs(wdset - rdset))

    passed = (passed .and. (ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm_world, mpi_err)
    endif

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI netCDF read 1D dataset', passed)
    endif

    deallocate(wdset)
    deallocate(rdset)
    deallocate(recvcounts)
    deallocate(n_size)

    call mpi_comm_finalise

end program test_mpi_netcdf_read_dataset_1d
