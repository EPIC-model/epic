! =============================================================================
!                       Test nearest algorithm
!
!   This unit test checks a = b - c across x MPI boundaries. The small parcels
!   'a' and 'c' should want to merge with the small parcel 'b'.
! =============================================================================
program test_mpi_nearest_8
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcels_mod, only : parcels
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels
    use parcel_nearest
    use mpi_environment
    use mpi_layout
    use mpi_timer
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

    call parcels%allocate(max_num_parcels)

    call parcel_setup

    parcels%local_num = n - 1
    parcels%total_num = 0

    call MPI_Allreduce(parcels%local_num,   &
                       parcels%total_num,   &
                       1,                   &
                       MPI_INTEGER,         &
                       MPI_SUM,             &
                       world%comm,          &
                       world%err)

    call find_nearest(parcels, isma, iclo, inva, n_merge, n_invalid)

    check_array(1) = parcels%local_num - n_invalid
    check_array(2) = n_merge

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, check_array, 2, MPI_INTEGER, MPI_SUM, world%root, world%comm, world%err)
    else
        call MPI_Reduce(check_array, check_array, 2, MPI_INTEGER, MPI_SUM, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        passed = (passed .and. (check_array(1) == parcels%total_num) .and. (check_array(2) == 200))

        call print_result_logical('Test MPI nearest algorithm: (4) a = b - c', passed)
    endif

    call nearest_win_deallocate

    call mpi_env_finalise


    contains

        subroutine cell_placement(l, i, j, k)
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

            ! small parcel c
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.24d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x
            parcels%position(2, l) = y - dx(2) * 0.42d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x
            parcels%position(2, l) = y + dx(2) * 0.46d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement

        subroutine parcel_setup
            integer :: i, j, k

            n = 1
            do k = box%lo(3)+1, box%lo(3)+1
                do j = box%lo(2), box%hi(2)
                    do i = box%lo(1), box%hi(1)
                        call cell_placement(n, i, j, k)
                    enddo
                enddo
            enddo

        end subroutine parcel_setup

end program test_mpi_nearest_8
