! =============================================================================
!                           Test parcel packing
!
!               This unit test checks the packing of parcels used
!               for sending the parcel around.
! =============================================================================
program test_mpi_parcel_pack
    use unit_test
    use constants
    use parcel_container
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

    n_parcels = 1000
    call parcel_alloc(n_parcels)

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
        parcels%B(5, n) = 12.0d0 + a
        parcels%volume(n) = 13.0d0 + a
        parcels%buoyancy(n) = 14.0d0 + a
#ifndef ENABLE_DRY_MODE
        parcels%humidity(n) = 15.0d0 + a
#endif
    enddo

    do m = 1, 100
        ! randomly select parcels to pack
        pid = -1
        do k = 1, n_pack
            call random_number(val)
            l = int(val * dble(n_parcels)) + 1
            i = findloc(pid, l, dim=1)
            do while (.not. i == 0)
                call random_number(val)
                l = int(val * dble(n_parcels)) + 1
                i = findloc(pid, l, dim=1)
            enddo
            pid(k) = l
        enddo

        allocate(buffer(n_par_attrib * n_pack))

        ! pack the parcels
        call parcel_pack(pid, n_pack, buffer)

        ! check the buffer of packed parcels
        do n = 1, n_pack

            l = pid(n)

            i = (n - 1) * n_par_attrib

            passed = (passed .and. (parcels%position(1, l) - buffer(i + IDX_X_POS) == zero))
            passed = (passed .and. (parcels%position(2, l) - buffer(i + IDX_Y_POS) == zero))
            passed = (passed .and. (parcels%position(3, l) - buffer(i + IDX_Z_POS) == zero))
            passed = (passed .and. (parcels%vorticity(1, l) - buffer(i + IDX_X_VOR) == zero))
            passed = (passed .and. (parcels%vorticity(2, l) - buffer(i + IDX_Y_VOR) == zero))
            passed = (passed .and. (parcels%vorticity(3, l) - buffer(i + IDX_Z_VOR) == zero))
            passed = (passed .and. (parcels%B(1, l) - buffer(i + IDX_B11) == zero))
            passed = (passed .and. (parcels%B(2, l) - buffer(i + IDX_B12) == zero))
            passed = (passed .and. (parcels%B(3, l) - buffer(i + IDX_B13) == zero))
            passed = (passed .and. (parcels%B(4, l) - buffer(i + IDX_B22) == zero))
            passed = (passed .and. (parcels%B(5, l) - buffer(i + IDX_B23) == zero))
            passed = (passed .and. (parcels%B(6, l) - buffer(i + IDX_B33) == zero))
            passed = (passed .and. (parcels%volume(l) - buffer(i + IDX_VOL) == zero))
            passed = (passed .and. (parcels%buoyancy(l) - buffer(i + IDX_BUO) == zero))
#ifndef ENABLE_DRY_MODE
            passed = (passed .and. (parcels%humidity(l) - buffer(i + IDX_HUM) == zero))
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
