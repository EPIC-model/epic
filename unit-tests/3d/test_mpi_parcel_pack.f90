! =============================================================================
!                           Test parcel packing
!
!               This unit test checks the packing of parcels used
!               for sending the parcel around.
! =============================================================================
program test_mpi_parcel_pack
    use unit_test
    use constants
    use parcels_mod, only : parcels
    use mpi_environment
    use mpi_timer
    implicit none

    logical :: passed = .true.
    integer :: n, m, k, i, l
    double precision :: a, val
    integer :: sk
    integer, parameter :: n_pack = 222
    integer, allocatable :: seed(:)
    integer :: pid(n_pack)
    double precision, allocatable:: buffer(:)

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = world%rank
    call random_seed(put=seed)

    parcels%local_num = 1000
    call parcels%allocate(parcels%local_num)

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

    do m = 1, 100
        ! randomly select parcels to pack
        pid = -1
        do k = 1, n_pack
            call random_number(val)
            l = int(val * dble(parcels%local_num)) + 1
            i = findloc(pid, l, dim=1)
            do while (.not. i == 0)
                call random_number(val)
                l = int(val * dble(parcels%local_num)) + 1
                i = findloc(pid, l, dim=1)
            enddo
            pid(k) = l
        enddo

        allocate(buffer(parcels%attr_num * n_pack))

        ! pack the parcels
        call parcels%pack(pid, n_pack, buffer)

        ! check the buffer of packed parcels
        do n = 1, n_pack

            l = pid(n)

            i = (n - 1) * parcels%attr_num

            passed = (passed .and. (parcels%position(1, l) - buffer(i + parcels%IDX_POS_BEG  ) == zero))
            passed = (passed .and. (parcels%position(2, l) - buffer(i + parcels%IDX_POS_BEG+1) == zero))
            passed = (passed .and. (parcels%position(3, l) - buffer(i + parcels%IDX_POS_BEG+2) == zero))
            passed = (passed .and. (parcels%vorticity(1, l) - buffer(i + parcels%IDX_VOR_BEG  ) == zero))
            passed = (passed .and. (parcels%vorticity(2, l) - buffer(i + parcels%IDX_VOR_BEG+1) == zero))
            passed = (passed .and. (parcels%vorticity(3, l) - buffer(i + parcels%IDX_VOR_BEG+2) == zero))
            passed = (passed .and. (parcels%B(1, l) - buffer(i + parcels%IDX_SHAPE_BEG  ) == zero))
            passed = (passed .and. (parcels%B(2, l) - buffer(i + parcels%IDX_SHAPE_BEG+1) == zero))
            passed = (passed .and. (parcels%B(3, l) - buffer(i + parcels%IDX_SHAPE_BEG+2) == zero))
            passed = (passed .and. (parcels%B(4, l) - buffer(i + parcels%IDX_SHAPE_BEG+3) == zero))
            passed = (passed .and. (parcels%B(5, l) - buffer(i + parcels%IDX_SHAPE_BEG+4) == zero))
            passed = (passed .and. (parcels%volume(l) - buffer(i + parcels%IDX_VOL) == zero))
            passed = (passed .and. (parcels%buoyancy(l) - buffer(i + parcels%IDX_BUO) == zero))
#ifndef ENABLE_DRY_MODE
            passed = (passed .and. (parcels%humidity(l) - buffer(i + parcels%IDX_HUM) == zero))
#endif
        enddo

        deallocate(buffer)
    enddo

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel pack', passed)
    endif

    call mpi_env_finalise

end program test_mpi_parcel_pack
