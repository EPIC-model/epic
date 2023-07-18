! =============================================================================
!                       Test nearest algorithm
!
!   This unit test checks a = b  C across MPI boundaries. The small
!   parcels 'a' and 'b' should want to merge.
! =============================================================================
program test_mpi_nearest_13
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels
    use parcel_nearest
    use mpi_communicator
    use mpi_layout
    use mpi_timer
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma, inva
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, check_array(2), n_invalid

    call mpi_comm_initialise

    passed = (passed .and. (world%err == 0))

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    parcel%lambda_max = five
    ! vmin = vcell / parcel%min_vratio
    parcel%min_vratio = 8.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call parcel_alloc(max_num_parcels)

    !
    ! Test 1: a=b C where 'a' and 'b' are on the same MPI rank; right diagonal MPI boundary
    !
    call parcel_setup(1)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(100)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (1) a = b  C', passed)
    endif

    !
    ! Test 2: a=b C where 'C' and 'b' are on the same MPI rank; right diagonal MPI boundary
    !
    call parcel_setup(2)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(100)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (2) a = b  C', passed)
    endif

    !
    ! Test 3: a=b C where 'a' and 'b' are on the same MPI rank; left diagonal MPI boundary
    !
    call parcel_setup(3)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(100)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (3) a = b  C', passed)
    endif

    !
    ! Test 4: a=b C where 'C' and 'b' are on the same MPI rank; left diagonal MPI boundary
    !
    call parcel_setup(4)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(100)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (4) a = b  C', passed)
    endif

    !
    ! Test 5: a=b C where 'a' and 'b' are on the same MPI rank
    !
    call parcel_setup(5)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (5) a = b  C', passed)
    endif

    !
    ! Test 6: a=b C where 'a' and 'b' are on the same MPI rank
    !
    call parcel_setup(6)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (6) a = b  C', passed)
    endif

    !
    ! Test 7: a=b C where 'C' and 'b' are on the same MPI rank
    !
    call parcel_setup(7)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (7) a = b  C', passed)
    endif

    !
    ! Test 8: a=b C where 'C' and 'b' are on the same MPI rank
    !
    call parcel_setup(8)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (8) a = b  C', passed)
    endif

    call nearest_win_deallocate

    call mpi_comm_finalise


    contains

        ! a = b  C where 'a' and 'b' are on the same MPI rank; right diagonal MPI boundary
        subroutine cell_placement_1(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x - dx(1) * 0.4d0
            parcels%position(2, l) = y - dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.36d0
            parcels%position(2, l) = y + dx(2) * 0.34d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_1

        ! a = b  C where 'C' and 'b' are on the same MPI rank; right diagonal MPI boundary
        subroutine cell_placement_2(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x - dx(1) * 0.42d0
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x + dx(1) * 0.28d0
            parcels%position(2, l) = y + dx(2) * 0.26d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_2

        ! a = b  C where 'a' and 'b' are on the same MPI rank; left diagonal MPI boundary
        subroutine cell_placement_3(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            ! small parcel b
            parcels%position(1, l) = x - dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x + dx(1) * 0.42d0
            parcels%position(2, l) = y - dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x - dx(1) * 0.36d0
            parcels%position(2, l) = y + dx(2) * 0.34d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_3

        ! a = b  C where 'C' and 'b' are on the same MPI rank; left diagonal MPI boundary
        subroutine cell_placement_4(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            ! small parcel b
            parcels%position(1, l) = x - dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y - dx(2) * 0.42d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x - dx(1) * 0.30d0
            parcels%position(2, l) = y + dx(2) * 0.28d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_4

        ! a = b  C where 'a' and 'b' are on the same MPI rank
        subroutine cell_placement_5(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            !
            ! in x-direction
            !

            ! small parcel b
            parcels%position(1, l) = x - dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x + dx(1) * 0.4d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x - dx(1) * 0.34d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            !
            ! in y-direction
            !

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.34d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_5

        ! a = b  C where 'a' and 'b' are on the same MPI rank
        subroutine cell_placement_6(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            !
            ! in x-direction
            !

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x - dx(1) * 0.4d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.34d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            !
            ! in y-direction
            !

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.34d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_6

        ! a = b  C where 'C' and 'b' are on the same MPI rank
        subroutine cell_placement_7(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            !
            ! in x-direction
            !

            ! small parcel b
            parcels%position(1, l) = x - dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x - dx(1) * 0.28d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            !
            ! in y-direction
            !

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.28d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_7

        ! a = b  C where 'C' and 'b' are on the same MPI rank
        subroutine cell_placement_8(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            !
            ! in x-direction
            !

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x + dx(1) * 0.28d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x - dx(1) * 0.44d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            !
            ! in y-direction
            !

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! big parcel C
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.28d0
            parcels%position(3, l) = z
            parcels%volume(l) = 1.1d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_8


        subroutine parcel_setup(num)
            integer, intent(in) :: num
            integer :: i, j, k

            n = 1
            do k = box%lo(3)+1, box%lo(3)+1
                do j = box%lo(2), box%hi(2)
                    do i = box%lo(1), box%hi(1)
                        select case (num)
                            case (1)
                                call cell_placement_1(n, i, j, k)
                            case (2)
                                call cell_placement_2(n, i, j, k)
                            case (3)
                                call cell_placement_3(n, i, j, k)
                            case (4)
                                call cell_placement_4(n, i, j, k)
                            case (5)
                                call cell_placement_5(n, i, j, k)
                            case (6)
                                call cell_placement_6(n, i, j, k)
                            case (7)
                                call cell_placement_7(n, i, j, k)
                            case (8)
                                call cell_placement_8(n, i, j, k)
                            case default
                                call mpi_exit_on_error("No valid parcel setup.")
                        end select
                    enddo
                enddo
            enddo

            n_parcels = n - 1
            n_total_parcels = 0

            call MPI_Allreduce(n_parcels,       &
                               n_total_parcels, &
                               1,               &
                               MPI_INTEGER,     &
                               MPI_SUM,         &
                               world%comm,      &
                               world%err)
        end subroutine parcel_setup

        subroutine check_result(n_true_merges)
            integer, intent(in) :: n_true_merges
            check_array(1) = n_parcels - n_invalid
            check_array(2) = n_merge

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE, check_array, 2, MPI_INTEGER, MPI_SUM, &
                                world%root, world%comm, world%err)
            else
                call MPI_Reduce(check_array, check_array, 2, MPI_INTEGER, MPI_SUM, &
                                world%root, world%comm, world%err)
            endif

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, &
                                world%root, world%comm, world%err)
            else
                call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, &
                                world%root, world%comm, world%err)
            endif

            if (world%rank == world%root) then
                passed = (passed .and. (check_array(1) == n_total_parcels) .and. (check_array(2) == n_true_merges))
            endif
        end subroutine check_result

end program test_mpi_nearest_13
