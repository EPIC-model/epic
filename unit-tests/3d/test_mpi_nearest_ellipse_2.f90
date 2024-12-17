! =============================================================================
!                       Test nearest algorithm
!
!   This unit test checks a = b - c across x MPI boundaries. The small parcels
!   'a' and 'c' should want to merge with the small parcel 'b'.
! =============================================================================
program test_mpi_nearest_8
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcels_mod, only : pc_type, parcels, bot_parcels, top_parcels
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, amin, max_num_parcels
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
    call register_timer('nearest MPI allreduce', nearest_allreduce_timer)

    parcel%lambda_max = five
    ! vmin = vcell / parcel%min_vratio
    parcel%min_vratio = 8.0d0
    parcel%min_aratio = 8.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call run_test(parcels, iz=box%lo(3)+1, rmin=vmin)
    bot_parcels%z_position = lower(3)
    call run_test(bot_parcels, iz=box%lo(3), rmin=amin)
    top_parcels%z_position = extent(3)
    call run_test(top_parcels, iz=box%lo(3), rmin=amin)

    call mpi_env_finalise


    contains

        subroutine run_test(pcont, iz, rmin)
            class(pc_type),   intent(inout) :: pcont
            integer,          intent(in)    :: iz
            double precision, intent(in)    :: rmin

            call pcont%allocate(max_num_parcels)

            call parcel_setup(pcont, iz, rmin)

            pcont%local_num = n - 1
            pcont%total_num = 0

            call MPI_Allreduce(pcont%local_num, &
                               pcont%total_num, &
                               1,               &
                               MPI_INTEGER,     &
                               MPI_SUM,         &
                               world%comm,      &
                               world%err)

            call find_nearest(pcont, isma, iclo, inva, n_merge, n_invalid)

            check_array(1) = pcont%local_num - n_invalid
            check_array(2) = n_merge

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE,   &
                                check_array,    &
                                2,              &
                                MPI_INTEGER,    &
                                MPI_SUM,        &
                                world%root,     &
                                world%comm,     &
                                world%err)
            else
                call MPI_Reduce(check_array,    &
                                check_array,    &
                                2,              &
                                MPI_INTEGER,    &
                                MPI_SUM,        &
                                world%root,     &
                                world%comm,     &
                                world%err)
            endif

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE,   &
                                passed,         &
                                1,              &
                                MPI_LOGICAL,    &
                                MPI_LAND,       &
                                world%root,     &
                                world%comm,     &
                                world%err)
            else
                call MPI_Reduce(passed,         &
                                passed,         &
                                1,              &
                                MPI_LOGICAL,    &
                                MPI_LAND,       &
                                world%root,     &
                                world%comm,     &
                                world%err)
            endif

            if (world%rank == world%root) then
                passed = (passed .and. (check_array(1) == pcont%total_num) .and. (check_array(2) == 200))

                call print_result_logical('Test MPI nearest algorithm: (4) a = b - c', passed)
            endif
        end subroutine run_test

        subroutine cell_placement(pcont, l, i, j, k, rmin)
            class(pc_type),   intent(inout) :: pcont
            integer,          intent(inout) :: l
            integer,          intent(in)    :: i, j, k
            double precision, intent(in)    :: rmin
            integer                         :: ix, iy, iz
            double precision                :: x, y, z

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            ! small parcel c
            pcont%position(1, l) = x
            pcont%position(2, l) = y - dx(2) * 0.24d0
            if (rmin /= amin) then
                pcont%position(3, l) = z
            endif
            pcont%volume(l) = 0.9d0 * rmin
            pcont%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel b
            pcont%position(1, l) = x
            pcont%position(2, l) = y - dx(2) * 0.42d0
            if (rmin /= amin) then
                pcont%position(3, l) = z
            endif
            pcont%volume(l) = 0.9d0 * rmin
            pcont%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel a
            pcont%position(1, l) = x
            pcont%position(2, l) = y + dx(2) * 0.46d0
            if (rmin /= amin) then
                pcont%position(3, l) = z
            endif
            pcont%volume(l) = 0.9d0 * rmin
            pcont%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement

        subroutine parcel_setup(pcont, iz, rmin)
            class(pc_type),   intent(inout) :: pcont
            integer,          intent(in)    :: iz
            double precision, intent(in)    :: rmin
            integer                         :: i, j

            n = 1
            do j = box%lo(2), box%hi(2)
                do i = box%lo(1), box%hi(1)
                    call cell_placement(pcont, n, i, j, iz, rmin)
                enddo
            enddo

        end subroutine parcel_setup

end program test_mpi_nearest_8
