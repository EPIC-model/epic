! =============================================================================
!                       Test nearest algorithm
!
!   This unit test checks a = b   c = d across MPI boundaries. The parcels
!   'a' and 'b' and the parcels 'c' and 'd' want to merge.
! =============================================================================
program test_mpi_nearest_18
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels
    use parcel_nearest
    use mpi_environment
    use mpi_layout
    use mpi_timer
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma, inva
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, check_array(2), n_invalid

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('nearest MPI barrier', nearest_barrier_timer)
    call register_timer('nearest MPI allreduce', nearest_allreduce_timer)

    parcel%lambda_max = five
    ! vmin = vcell / parcel%min_vratio
    parcel%min_vratio = 8.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call parcel_alloc(max_num_parcels)

    !
    ! Test 1: a = b  c = d; corner setup
    !
    call parcel_setup(1)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (1) a = b  c = d', passed)
    endif

    !
    ! Test 2: a = b  c = d; corner setup (90 degree rotation of test 1)
    !
    call parcel_setup(2)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(200)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (2) a = b  c = d', passed)
    endif

    !
    ! Test 3: a = b  c = d with a = b across MPI boundary
    !
    call parcel_setup(3)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(400)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (3) a = b  c = d', passed)
    endif

    !
    ! Test 4: a = b  c = d with c = d across MPI boundary
    !
    call parcel_setup(4)

    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    call check_result(400)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm: (4) a = b  c = d', passed)
    endif

    call nearest_win_deallocate

    call mpi_env_finalise


    contains

        ! a = b  c = d; corner setup
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

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.32d0
            parcels%position(2, l) = y - dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.45d0
            parcels%position(2, l) = y - dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y + dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x - dx(1) * 0.32d0
            parcels%position(2, l) = y + dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

        end subroutine cell_placement_1

        ! a = b  c = d; corner setup (90 degree rotation of test 1)
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

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.45d0
            parcels%position(2, l) = y - dx(2) * 0.3d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.45d0
            parcels%position(2, l) = y - dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y + dx(2) * 0.3d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y + dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

        end subroutine cell_placement_2

        ! a = b  c = d with a = b across MPI boundary
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

            !
            ! in x-direction
            !

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.45d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.3d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x - dx(1) * 0.2d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1


            !
            ! in y-direction
            !

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.45d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.45d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.3d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.2d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

        end subroutine cell_placement_3

        ! a = b  c = d with c = d across MPI boundary
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

            !
            ! in x-direction
            !

            ! small parcel a
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.45d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x + dx(1) * 0.3d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x + dx(1) * 0.2d0
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1


            !
            ! in y-direction
            !

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.45d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.45d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.3d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.2d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%theta(l) = l + world%rank * 100

            l = l + 1

        end subroutine cell_placement_4

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

end program test_mpi_nearest_18
