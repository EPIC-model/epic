! =============================================================================
!                           Test parcel deletion
!
!               This unit test checks the deletion of parcels.
! =============================================================================
program test_mpi_parcel_delete
    use unit_test
    use constants
    use parcels_mod, only : parcels
    use mpi_environment
    use mpi_timer
    use merge_sort
    implicit none

    logical :: passed = .true.
    integer :: n, m, k, i, l
    double precision :: a, val
    integer :: sk
    integer, parameter :: n_del = 222
    integer, allocatable :: seed(:)
    integer :: invalid(0:n_del)
    integer :: ii(1000 - n_del)


    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = world%rank
    call random_seed(put=seed)

    do m = 1, 100
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

        ! randomly select parcels to delete
        invalid = -1
        do k = 1, n_del
            call random_number(val)
            l = int(val * dble(parcels%local_num)) + 1
            i = findloc(invalid, l, dim=1)
            do while (.not. i == 0)
                call random_number(val)
                l = int(val * dble(parcels%local_num)) + 1
                i = findloc(invalid, l, dim=1)
            enddo
            invalid(k) = l
        enddo

        call msort(invalid(1:n_del), ii(1:n_del))

        ! delete the parcels
        call parcels%delete(invalid, n_del)

        ! sort as the delete algorithm fills gaps with parcels from the back
        call msort(parcels%position(1, 1:parcels%local_num), ii)

        ! check remaining parcels
        k = 1
        i = 1
        do n = 1, parcels%local_num
            do while (i == invalid(k))
                i = i + 1
                k = k + 1
                if (k > n_del) then
                    k = 1
                    i = -1
                endif
            enddo

            if (i > 0) then
                a = dble(i-1) * 100.0d0

                ! Note: x-posiiton is sorted that is why we access with n and not ii(n)
                passed = (passed .and. (parcels%position(1, n) - (1.0d0 + a) == zero))
                passed = (passed .and. (parcels%position(2, ii(n)) - (2.0d0 + a) == zero))
                passed = (passed .and. (parcels%position(3, ii(n)) - (3.0d0 + a) == zero))
                passed = (passed .and. (parcels%vorticity(1, ii(n)) - (4.0d0 + a) == zero))
                passed = (passed .and. (parcels%vorticity(2, ii(n)) - (5.0d0 + a) == zero))
                passed = (passed .and. (parcels%vorticity(3, ii(n)) - (6.0d0 + a) == zero))
                passed = (passed .and. (parcels%B(1, ii(n)) - (7.0d0 + a) == zero))
                passed = (passed .and. (parcels%B(2, ii(n)) - (8.0d0 + a) == zero))
                passed = (passed .and. (parcels%B(3, ii(n)) - (9.0d0 + a) == zero))
                passed = (passed .and. (parcels%B(4, ii(n)) - (10.0d0 + a) == zero))
                passed = (passed .and. (parcels%B(5, ii(n)) - (11.0d0 + a) == zero))
                passed = (passed .and. (parcels%volume(ii(n)) - (12.0d0 + a) == zero))
                passed = (passed .and. (parcels%buoyancy(ii(n)) - (13.0d0 + a) == zero))
#ifndef ENABLE_DRY_MODE
                passed = (passed .and. (parcels%humidity(ii(n)) - (14.0d0 + a) == zero))
#endif
                i = i + 1
            endif
        enddo

        call parcels%deallocate
    enddo

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel delete', passed)
    endif

    call mpi_env_finalise

end program test_mpi_parcel_delete
