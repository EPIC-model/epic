! =============================================================================
!                           Test parcel unpacking
!
!               This unit test checks the unpacking of parcels used
!               when receiving parcels.
! =============================================================================
program test_mpi_parcel_unpack
    use unit_test
    use constants
    use parcel_container
    use mpi_communicator
    use mpi_timer
    implicit none

    logical :: passed = .true.
    integer :: n, i
    double precision :: a
    integer, parameter :: n_unpack = 222
    double precision, allocatable:: buffer(:)

    call mpi_comm_initialise

    passed = (passed .and. (comm%err == 0))

    n_parcels = 1000 - n_unpack
    call parcel_alloc(n_parcels + n_unpack)

    ! initialise parcels
    do n = 1, n_parcels
        a = dble(n-1) * 100.0d0
        parcels%position(1, n) = 1.0d0 + a
        parcels%position(2, n) = 2.0d0 + a
        parcels%position(3, n) = 3.0d0 + a
        parcels%vorticity(1, n) = 4.0d0 + a
        parcels%vorticity(2, n) = 5.0d0 + a
        parcels%vorticity(3, n) = 6.0d0 + a
        parcels%B(1, n) = 7.0d0 + a
        parcels%B(2, n) = 8.0d0 + a
        parcels%B(3, n) = 9.0d0 + a
        parcels%B(4, n) = 10.0d0 + a
        parcels%B(5, n) = 11.0d0 + a
        parcels%volume(n) = 12.0d0 + a
        parcels%buoyancy(n) = 13.0d0 + a
#ifndef ENABLE_DRY_MODE
        parcels%humidity(n) = 14.0d0 + a
#endif
    enddo

    allocate(buffer(n_par_attrib * n_unpack))

    do n = 1, n_unpack
        a = dble(n+n_parcels-1) * 100.0d0

        i = (n - 1) * n_par_attrib

        buffer(i + IDX_X_POS) = 1.0d0 + a
        buffer(i + IDX_Y_POS) = 2.0d0 + a
        buffer(i + IDX_Z_POS) = 3.0d0 + a
        buffer(i + IDX_X_VOR) = 4.0d0 + a
        buffer(i + IDX_Y_VOR) = 5.0d0 + a
        buffer(i + IDX_Z_VOR) = 6.0d0 + a
        buffer(i + IDX_B11) = 7.0d0 + a
        buffer(i + IDX_B12) = 8.0d0 + a
        buffer(i + IDX_B13) = 9.0d0 + a
        buffer(i + IDX_B22) = 10.0d0 + a
        buffer(i + IDX_B23) = 11.0d0 + a
        buffer(i + IDX_VOL) = 12.0d0 + a
        buffer(i + IDX_BUO) = 13.0d0 + a
#ifndef ENABLE_DRY_MODE
        buffer(i + IDX_HUM) = 14.0d0 + a
#endif
    enddo

    ! pack the parcels
    call parcel_unpack(n_unpack, buffer)

    passed = (passed .and. (n_parcels == 1000))

    ! check the parcel container
    do n = 1, n_parcels
        a = dble(n-1) * 100.0d0

        passed = (passed .and. (parcels%position(1, n) == 1.0d0 + a))
        passed = (passed .and. (parcels%position(2, n) == 2.0d0 + a))
        passed = (passed .and. (parcels%position(3, n) == 3.0d0 + a))
        passed = (passed .and. (parcels%vorticity(1, n) == 4.0d0 + a))
        passed = (passed .and. (parcels%vorticity(2, n) == 5.0d0 + a))
        passed = (passed .and. (parcels%vorticity(3, n) == 6.0d0 + a))
        passed = (passed .and. (parcels%B(1, n) == 7.0d0 + a))
        passed = (passed .and. (parcels%B(2, n) == 8.0d0 + a))
        passed = (passed .and. (parcels%B(3, n) == 9.0d0 + a))
        passed = (passed .and. (parcels%B(4, n) == 10.0d0 + a))
        passed = (passed .and. (parcels%B(5, n) == 11.0d0 + a))
        passed = (passed .and. (parcels%volume(n) == 12.0d0 + a))
        passed = (passed .and. (parcels%buoyancy(n) == 13.0d0 + a))
#ifndef ENABLE_DRY_MODE
        passed = (passed .and. (parcels%humidity(n) == 14.0d0 + a))
#endif
    enddo

    deallocate(buffer)

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI parcel unpack', passed)
    endif

    call mpi_comm_finalise

end program test_mpi_parcel_unpack
