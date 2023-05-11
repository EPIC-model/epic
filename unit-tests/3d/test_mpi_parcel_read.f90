! =============================================================================
!                           Test parallel parcel reading
!
!       This unit test checks the reading of parcels in parallel.
! =============================================================================
program test_mpi_parcel_read
    use unit_test
    use constants, only : zero
    use parcel_container
    use parameters, only : set_max_num_parcels
    use parcel_netcdf
    use mpi_communicator
    use mpi_timer
    implicit none

    logical              :: passed = .true.
    double precision     :: res
    integer              :: n, start_index, n_parcels_before

    call set_max_num_parcels(100000000)

    call mpi_comm_initialise

    passed = (passed .and. (comm%err == 0))

    call register_timer('parcel I/O', parcel_io_timer)

    !
    ! write parcels first
    !

    n_parcels = 10 + (comm%rank + 1)
    n_parcels_before = n_parcels
    n_total_parcels = 11 * comm%size + (comm%size * (comm%size - 1)) / 2

    call parcel_alloc(n_total_parcels)

    ! fill with 1 to n_total_parcels
    start_index = 0
    do n = 0, comm%rank-1
        start_index = start_index + 10 + (n + 1)
    enddo

    do n = 1, n_parcels
        parcels%position(:, n) = start_index + n
        parcels%B(:, n) = start_index + n
        parcels%volume(n) = start_index + n
        parcels%vorticity(:, n) = start_index + n
        parcels%theta(n) = start_index + n

#ifndef ENABLE_DRY_MODE
        parcels%qv(n) = start_index + n
        parcels%ql(n) = start_index + n
#endif
    enddo

    call create_netcdf_parcel_file('nctest', .true., .false.)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_parcels(t = 10.0d0)

    passed = (passed .and. (ncerr == 0))

    !
    ! read parcels now
    !

    parcels%position = 0
    parcels%B = 0
    parcels%volume = 0
    parcels%vorticity = 0
    parcels%theta = 0
#ifndef ENABLE_DRY_MODE
    parcels%qv = 0
    parcels%ql = 0
#endif

    call read_netcdf_parcels('nctest_0000000001_parcels.nc')

    passed = (passed .and. (n_parcels == n_parcels_before))

    if (passed) then
        do n = 1, n_parcels
            res = dble(start_index + n)
            passed = (passed .and. (abs(parcels%position(1, n) - res) == zero))
            passed = (passed .and. (abs(parcels%position(2, n) - res) == zero))
            passed = (passed .and. (abs(parcels%position(3, n) - res) == zero))
            passed = (passed .and. (abs(parcels%B(1, n) - res) == zero))
            passed = (passed .and. (abs(parcels%B(2, n) - res) == zero))
            passed = (passed .and. (abs(parcels%B(3, n) - res) == zero))
            passed = (passed .and. (abs(parcels%B(4, n) - res) == zero))
            passed = (passed .and. (abs(parcels%B(5, n) - res) == zero))
            passed = (passed .and. (abs(parcels%volume(n) - res) == zero))
            passed = (passed .and. (abs(parcels%vorticity(1, n) - res) == zero))
            passed = (passed .and. (abs(parcels%vorticity(2, n) - res) == zero))
            passed = (passed .and. (abs(parcels%vorticity(3, n) - res) == zero))
            passed = (passed .and. (abs(parcels%theta(n) - res) == zero))
#ifndef ENABLE_DRY_MODE
            passed = (passed .and. (abs(parcels%qv(n) - res) == zero))
            passed = (passed .and. (abs(parcels%ql(n) - res) == zero))
#endif
        enddo
    endif

    call parcel_dealloc

    call delete_netcdf_file(ncfname='nctest_0000000001_parcels.nc')

    passed = (passed .and. (ncerr == 0))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI parcel read', passed)
    endif

    call mpi_comm_finalise

end program test_mpi_parcel_read
