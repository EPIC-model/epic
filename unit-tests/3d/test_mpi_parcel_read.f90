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
    use netcdf_utils
    use mpi_environment
    use mpi_timer
    implicit none

    logical              :: passed = .true.
    double precision     :: res
    integer              :: n, start_index, n_parcels_before 

    call set_max_num_parcels(100000000)

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    call register_timer('parcel I/O', parcel_io_timer)

    !
    ! write parcels first
    !

    n_parcels = 10 + (world%rank + 1)
    n_parcels_before = n_parcels
    n_total_parcels = 11 * world%size + (world%size * (world%size - 1)) / 2

    call parcel_alloc(n_total_parcels)

    ! fill with 1 to n_total_parcels
    start_index = 0
    do n = 0, world%rank-1
        start_index = start_index + 10 + (n + 1)
    enddo

    do n = 1, n_parcels
        parcels%position(:, n) = start_index + n
        parcels%B(:, n) = start_index + n
        parcels%volume(n) = start_index + n
        parcels%vorticity(:, n) = start_index + n
        parcels%buoyancy(n) = start_index + n

#ifndef ENABLE_DRY_MODE
        parcels%humidity(n) = start_index + n
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
    parcels%buoyancy = 0
#ifndef ENABLE_DRY_MODE
    parcels%humidity = 0
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
            passed = (passed .and. (abs(parcels%buoyancy(n) - res) == zero))
#ifndef ENABLE_DRY_MODE
            passed = (passed .and. (abs(parcels%humidity(n) - res) == zero))
#endif
        enddo
    endif

    call parcel_dealloc

    call delete_netcdf_file(ncfname='nctest_0000000001_parcels.nc')

    passed = (passed .and. (ncerr == 0))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel read', passed)
    endif

    call mpi_env_finalise

end program test_mpi_parcel_read
