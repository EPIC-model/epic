! =============================================================================
!                           Test parallel parcel writing
!
!       This unit test checks the writing of parcels in parallel.
! =============================================================================
program test_mpi_parcel_write
    use unit_test
    use constants, only : zero
    use parcel_container
    use parcel_netcdf
    use mpi_communicator
    use timer
    implicit none

    logical :: passed = .true.

    call mpi_comm_initialise

    passed = (passed .and. (mpi_err == 0))

    call register_timer('parcel I/O', parcel_io_timer)

    n_parcels = 10 + (mpi_rank + 1)
    n_total_parcels = 11 * mpi_size + (mpi_size * (mpi_size - 1)) / 2

    call parcel_alloc(n_parcels)

    parcels%position(:, 1:n_parcels) = mpi_rank
    parcels%B(:, 1:n_parcels) = mpi_rank
    parcels%volume(1:n_parcels) = mpi_rank
    parcels%vorticity(:, 1:n_parcels) = mpi_rank
    parcels%buoyancy(1:n_parcels) = mpi_rank

#ifndef ENABLE_DRY_MODE
    parcels%humidity(1:n_parcels) = mpi_rank
#endif

    call create_netcdf_parcel_file('nctest', .true., .false.)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_parcels(t = 10.0d0)

    passed = (passed .and. (ncerr == 0))

    call parcel_dealloc

    call delete_netcdf_file(ncfname='nctest_0000000001_parcels.nc')

    passed = (passed .and. (ncerr == 0))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm_world, mpi_err)
    endif

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI parcel write', passed)
    endif

    call mpi_comm_finalise

end program test_mpi_parcel_write
