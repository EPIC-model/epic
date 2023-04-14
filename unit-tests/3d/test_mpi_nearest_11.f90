! =============================================================================
!                       Test nearest algorithm
!
!   This unit test checks a = b - c across diagonal MPI boundaries. The small
!   parcels 'a' and 'c' should want to merge with the small parcel 'b'.
! =============================================================================
program test_mpi_nearest_11
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels
    use parcel_nearest
    use mpi_communicator
    use mpi_layout
    use mpi_timer
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma, inva
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, check_array(2), n_invalid

    call mpi_comm_initialise

    passed = (passed .and. (comm%err == 0))

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

    call update_parameters

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call parcel_alloc(max_num_parcels)

    call parcel_setup

    n_parcels = n - 1
    n_total_parcels = 0

    call MPI_Allreduce(n_parcels,       &
                       n_total_parcels, &
                       1,               &
                       MPI_INTEGER,     &
                       MPI_SUM,         &
                       comm%world,      &
                       comm%err)


    call find_nearest(isma, iclo, inva, n_merge, n_invalid)

    check_array(1) = n_parcels - n_invalid
    check_array(2) = n_merge

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, check_array, 2, MPI_INTEGER, MPI_SUM, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(check_array, check_array, 2, MPI_INTEGER, MPI_SUM, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        passed = (passed .and. (check_array(1) == n_total_parcels) .and. (check_array(2) == 200))

        call print_result_logical('Test MPI nearest algorithm: (7) a = b - c', passed)
    endif

    call mpi_comm_finalise


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

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + comm%rank * 100
            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.4d0
            parcels%position(2, l) = y - dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + comm%rank * 100
            l = l + 1

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.35d0
            parcels%position(2, l) = y + dx(2) * 0.33d0
            parcels%position(3, l) = z
            parcels%volume(l) = 0.9d0 * vmin
            parcels%buoyancy(l) = l + comm%rank * 100
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

end program test_mpi_nearest_11
