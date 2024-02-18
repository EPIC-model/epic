! =============================================================================
!                       Test nearest algorithm
!
!    This unit test checks dual links (a = b) across MPI boundaries.
! =============================================================================
program test_mpi_nearest_ellipse_1
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcels_mod, only : bot_parcels
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, amin, max_num_parcels
    use parcel_nearest
    use mpi_environment
    use mpi_layout
    use mpi_timer
    use mpi_utils, only : mpi_stop
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma, inva
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, check_array(2), n_invalid
    integer :: r, is, ic

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    parcel%lambda_max = five
    ! amin = acell / parcel%min_aratio
    parcel%min_aratio = 8.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call bot_parcels%allocate(max_num_parcels)

    call parcel_setup

    bot_parcels%local_num = n - 1
    bot_parcels%total_num = 0

    call MPI_Allreduce(bot_parcels%local_num,   &
                       bot_parcels%total_num,   &
                       1,                       &
                       MPI_INTEGER,             &
                       MPI_SUM,                 &
                       world%comm,              &
                       world%err)


    call find_nearest(bot_parcels, isma, iclo, inva, n_merge, n_invalid)

    ! To print out the result enable the following lines:
    do r = 0, world%size-1
        if (r == world%rank) then
            do n = 1, n_merge
                is = isma(n)
                ic = iclo(n)
                print *, world%rank, bot_parcels%position(1, is),                                &
                                     bot_parcels%position(2, is), int(bot_parcels%buoyancy(is)), &
                                     bot_parcels%position(1, ic),                                &
                                     bot_parcels%position(2, ic), int(bot_parcels%buoyancy(ic))
            enddo
        endif
        call MPI_Barrier(world%comm, world%err)
    enddo
    call mpi_stop

    check_array(1) = bot_parcels%local_num - n_invalid
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
        passed = (passed .and. (check_array(1) == bot_parcels%total_num) .and. (check_array(2) == 200))

        call print_result_logical('Test MPI nearest ellipse algorithm: a = b', passed)
    endif

    call nearest_win_deallocate

    call mpi_env_finalise


    contains

        subroutine cell_placement(l, i, j)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j
            integer                :: ix, iy, m
            double precision       :: x, y, z

            ix = i
            iy = j

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3)

            ! parcel in centre
            bot_parcels%position(1, l) = x
            bot_parcels%position(2, l) = y
            bot_parcels%z_position = z
            bot_parcels%area(l) = 1.1d0 * amin
            bot_parcels%buoyancy(l) = l + world%rank * 100

            l = l + 1

            do m = -1, 1, 2
                bot_parcels%position(1, l) = x + dble(m) * dx(1) * 0.45
                bot_parcels%position(2, l) = y
                bot_parcels%z_position = z
                bot_parcels%area(l) = 0.9d0 * amin
                bot_parcels%buoyancy(l) = l + world%rank * 100
                l = l + 1
            enddo

            do m = -1, 1, 2
                bot_parcels%position(1, l) = x
                bot_parcels%position(2, l) = y + dble(m) * dx(2) * 0.45
                bot_parcels%z_position = z
                bot_parcels%area(l) = 0.9d0 * amin
                bot_parcels%buoyancy(l) = l + world%rank * 100
                l = l + 1
            enddo

        end subroutine cell_placement

        subroutine parcel_setup
!             use mpi_utils, only : mpi_stop
            integer :: i, j!, r

            n = 1
            do j = box%lo(2), box%hi(2)
                do i = box%lo(1), box%hi(1)
                    call cell_placement(n, i, j)
                enddo
            enddo

!                         ! To print out the parcel setups enable the following lines:
!             do r = 0, world%size-1
!                 if (r == world%rank) then
!                     do i = 1, n - 1
!                         print *, bot_parcels%position(1, i), bot_parcels%position(2, i), i + 100 * world%rank, world%rank
!                     enddo
!                 endif
!                 call MPI_Barrier(world%comm, world%err)
!             enddo
!             call MPI_Barrier(world%comm, world%err)
!             call mpi_stop

        end subroutine parcel_setup

end program test_mpi_nearest_ellipse_1