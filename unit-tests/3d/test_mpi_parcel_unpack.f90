! =============================================================================
!                           Test parcel unpacking
!
!               This unit test checks the unpacking of parcels used
!               when receiving parcels.
! =============================================================================
program test_mpi_parcel_unpack
    use unit_test
    use constants
    use parcels_mod, only : parcels
    use mpi_environment
    use mpi_timer
    implicit none

    logical :: passed = .true.
    integer :: n, i
    double precision :: a
    integer, parameter :: n_unpack = 222
    double precision, allocatable:: buffer(:)

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    parcels%local_num = 1000 - n_unpack
    call parcels%allocate(parcels%local_num + n_unpack)

    ! initialise parcels
    do n = 1, parcels%local_num
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

    allocate(buffer(parcels%attr_num * n_unpack))

    do n = 1, n_unpack
        a = dble(n+parcels%local_num-1) * 100.0d0

        i = (n - 1) * parcels%attr_num

        buffer(i + parcels%IDX_POS_BEG  ) = 1.0d0 + a
        buffer(i + parcels%IDX_POS_BEG+1) = 2.0d0 + a
        buffer(i + parcels%IDX_POS_BEG+2) = 3.0d0 + a
        buffer(i + parcels%IDX_VOR_BEG  ) = 4.0d0 + a
        buffer(i + parcels%IDX_VOR_BEG+1) = 5.0d0 + a
        buffer(i + parcels%IDX_VOR_BEG+2) = 6.0d0 + a
        buffer(i + parcels%IDX_SHAPE_BEG  ) = 7.0d0 + a
        buffer(i + parcels%IDX_SHAPE_BEG+1) = 8.0d0 + a
        buffer(i + parcels%IDX_SHAPE_BEG+2) = 9.0d0 + a
        buffer(i + parcels%IDX_SHAPE_BEG+3) = 10.0d0 + a
        buffer(i + parcels%IDX_SHAPE_BEG+4) = 11.0d0 + a
        buffer(i + parcels%IDX_VOL) = 12.0d0 + a
        buffer(i + parcels%IDX_BUO) = 13.0d0 + a
#ifndef ENABLE_DRY_MODE
        buffer(i + parcels%IDX_HUM) = 14.0d0 + a
#endif
    enddo

    ! pack the parcels
    call parcels%unpack(n_unpack, buffer)

    passed = (passed .and. (parcels%local_num == 1000))

    ! check the parcel container
    do n = 1, parcels%local_num
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

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel unpack', passed)
    endif

    call mpi_env_finalise

end program test_mpi_parcel_unpack
