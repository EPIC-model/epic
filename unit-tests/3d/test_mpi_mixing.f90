! =============================================================================
!                       Test mixing algorithm
! =============================================================================
program test_mpi_mixing
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcel_container
    use parcels_mod
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels, amin
    use parcel_nearest
    use mpi_environment
    use mpi_layout
    use mpi_timer
    use parcel_mixing, only : mix_parcels, mixing_timer
    use mpi_utils, only : mpi_exit_on_error
    use mpi_utils, only : mpi_stop
    implicit none

    logical :: passed = .true.
    integer :: l

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call register_timer('parcel mixing', mixing_timer)

    parcel%lambda_max = five
    ! vmin = vcell / parcel%min_vratio
    parcel%min_vratio = 8.0d0
    parcel%min_aratio = 8.0d0

    parcel%n_surf_per_cell = 4

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call parcels%allocate(max_num_parcels)
    call top_parcels%allocate(nx * ny * 10)
    call bot_parcels%allocate(nx * ny * 10)

    do l = 1, 4
        call parcel_setup(l)

        call mix_parcels

        call check_result(l)

    enddo

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel mixing', passed)
    endif

    call mpi_env_finalise

    contains

        subroutine cell_placement_interior_1(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            double precision       :: x, y, z

            x = lower(1) + (0.5d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            parcels%position(1, l) = x
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.5d0 * vmin
            parcels%buoyancy(l) = 2.0d0
            l = l + 1

        end subroutine cell_placement_interior_1

        subroutine cell_placement_surface_1(psurf, l, i, j, k)
            type(ellipse_pc_type), intent(inout) :: psurf
            integer,               intent(inout) :: l
            integer,               intent(in)    :: i, j, k
            double precision                     :: x, y, z

            x = lower(1) + (0.5d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            psurf%position(1, l) = x
            psurf%position(2, l) = y
            psurf%z_position = z
            psurf%volume(l) = 0.5d0 * vmin
            psurf%area(l) = 0.5d0 * amin
            psurf%buoyancy(l) = 1.0d0
            l = l + 1

        end subroutine cell_placement_surface_1

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine cell_placement_interior_2(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            double precision       :: x, y, z

            x = lower(1) + (0.9d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            parcels%position(1, l) = x
            parcels%position(2, l) = y
            parcels%position(3, l) = z
            parcels%volume(l) = 1.5d0 * vmin
            parcels%buoyancy(l) = 2.0d0
            l = l + 1

        end subroutine cell_placement_interior_2

        subroutine cell_placement_surface_2(psurf, l, i, j, k)
            type(ellipse_pc_type), intent(inout) :: psurf
            integer,               intent(inout) :: l
            integer,               intent(in)    :: i, j, k
            double precision                     :: x, y, z

            x = lower(1) + (0.1d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            psurf%position(1, l) = x
            psurf%position(2, l) = y
            psurf%z_position = z
            psurf%volume(l) = 0.5d0 * vmin
            psurf%area(l) = 0.5d0 * amin
            psurf%buoyancy(l) = 1.0d0
            l = l + 1

        end subroutine cell_placement_surface_2

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine cell_placement_interior_3(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            double precision       :: x, y, z
            double precision       :: dd

            call cell_placement_interior_1(l, i, j, k)

            x = lower(1) + (0.5d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            if (k == 0) then
                dd = dx(3) * 0.001d0
            else if (k == nz) then
                dd = -dx(3) * 0.001d0
            else
                dd = 0.0d0
            endif

            parcels%position(1, l) = x
            parcels%position(2, l) = y
            parcels%position(3, l) = z + dd
            parcels%volume(l) = 0.5d0 * vmin
            parcels%buoyancy(l) = 1.0d0
            l = l + 1
        end subroutine cell_placement_interior_3

        subroutine cell_placement_surface_3(psurf, l, i, j, k)
            type(ellipse_pc_type), intent(inout) :: psurf
            integer,               intent(inout) :: l
            integer,               intent(in)    :: i, j, k
            double precision                     :: x, y, z

            call cell_placement_surface_1(psurf, l, i, j, k)

            x = lower(1) + (0.5d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            psurf%position(1, l) = x
            psurf%position(2, l) = y
            psurf%z_position = z
            psurf%volume(l) = 1.5d0 * vmin
            psurf%area(l) = 1.5d0 * amin
            psurf%buoyancy(l) = 2.0d0
            l = l + 1

        end subroutine cell_placement_surface_3

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine cell_placement_interior_4(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            double precision       :: x, y, z
            double precision       :: dd

            call cell_placement_interior_2(l, i, j, k)

            x = lower(1) + (0.9d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            if (k == 0) then
                dd = dx(3) * 0.001d0
            else if (k == nz) then
                dd = -dx(3) * 0.001d0
            else
                dd = 0.0d0
            endif

            parcels%position(1, l) = x
            parcels%position(2, l) = y
            parcels%position(3, l) = z + dd
            parcels%volume(l) = 0.5d0 * vmin
            parcels%buoyancy(l) = 1.0d0
            l = l + 1

        end subroutine cell_placement_interior_4

        subroutine cell_placement_surface_4(psurf, l, i, j, k)
            type(ellipse_pc_type), intent(inout) :: psurf
            integer,               intent(inout) :: l
            integer,               intent(in)    :: i, j, k
            double precision                     :: x, y, z

            call cell_placement_surface_2(psurf, l, i, j, k)

            x = lower(1) + (0.1d0 + dble(i)) * dx(1)
            y = lower(2) + (0.5d0 + dble(j)) * dx(2)
            z = lower(3) + dble(k) * dx(3)

            psurf%position(1, l) = x
            psurf%position(2, l) = y
            psurf%z_position = z
            psurf%volume(l) = 1.5d0 * vmin
            psurf%area(l) = 1.5d0 * amin
            psurf%buoyancy(l) = 2.0d0
            l = l + 1
        end subroutine cell_placement_surface_4


        subroutine parcel_setup(num)
            integer, intent(in) :: num
            integer :: i, j, np, nt, nb

            np = 1
            nb = 1
            nt = 1
            do j = box%lo(2), box%hi(2)
                do i = box%lo(1), box%hi(1)
                    select case (num)
                        case (1)
                            call cell_placement_interior_1(np, i, j, 0)
                            call cell_placement_interior_1(np, i, j, nz)

                            call cell_placement_surface_1(bot_parcels, nb, i, j, 0)
                            call cell_placement_surface_1(top_parcels, nt, i, j, nz)

                        case (2)
                            call cell_placement_interior_2(np, i, j, 0)
                            call cell_placement_interior_2(np, i, j, nz)

                            call cell_placement_surface_2(bot_parcels, nb, i, j, 0)
                            call cell_placement_surface_2(top_parcels, nt, i, j, nz)

                        case (3)
                            call cell_placement_interior_3(np, i, j, 0)
                            call cell_placement_interior_3(np, i, j, nz)

                            call cell_placement_surface_3(bot_parcels, nb, i, j, 0)
                            call cell_placement_surface_3(top_parcels, nt, i, j, nz)

                        case (4)
                            call cell_placement_interior_4(np, i, j, 0)
                            call cell_placement_interior_4(np, i, j, nz)

                            call cell_placement_surface_4(bot_parcels, nb, i, j, 0)
                            call cell_placement_surface_4(top_parcels, nt, i, j, nz)

                        case default
                            call mpi_exit_on_error("No valid parcel setup.")
                    end select
                enddo
            enddo

            parcels%local_num = np - 1
            parcels%total_num = 0

            call MPI_Allreduce(parcels%local_num,   &
                               parcels%total_num,   &
                               1,                   &
                               MPI_INTEGER,         &
                               MPI_SUM,             &
                               world%comm,          &
                               world%err)

            bot_parcels%local_num = nb - 1
            bot_parcels%total_num = 0

            call MPI_Allreduce(bot_parcels%local_num,   &
                               bot_parcels%total_num,   &
                               1,                       &
                               MPI_INTEGER,             &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            top_parcels%local_num = nt - 1
            top_parcels%total_num = 0

            call MPI_Allreduce(top_parcels%local_num,   &
                               top_parcels%total_num,   &
                               1,                       &
                               MPI_INTEGER,             &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

        end subroutine parcel_setup

        subroutine check_result(num)
            integer,         intent(in) :: num
            double precision            :: local_sum
            double precision            :: global_sum(3)

            local_sum = sum(parcels%buoyancy(1:parcels%local_num))
            call MPI_Allreduce(local_sum,               &
                               global_sum(1),           &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            local_sum = sum(top_parcels%buoyancy(1:top_parcels%local_num))
            call MPI_Allreduce(local_sum,               &
                               global_sum(2),           &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            local_sum = sum(bot_parcels%buoyancy(1:bot_parcels%local_num))
            call MPI_Allreduce(local_sum,               &
                               global_sum(3),           &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            select case (num)
                case (1)
                    ! buoyancy per parcel = (1.5 * 2.0 + 0.5 * 1.0)
                    passed = (passed .and. (global_sum(1) == 350.0d0))
                    passed = (passed .and. (global_sum(2) == 175.0d0))
                    passed = (passed .and. (global_sum(3) == 175.0d0))

                case (2)
                    ! buoyancy per parcel = (1.5 * 2.0 + 0.5 * 1.0)
                    passed = (passed .and. (global_sum(1) == 350.0d0))
                    passed = (passed .and. (global_sum(2) == 175.0d0))
                    passed = (passed .and. (global_sum(3) == 175.0d0))

                case (3)
                    ! buoyancy per parcel = (1.5 * 2.0 + 0.5 * 1.0)
                    passed = (passed .and. (global_sum(1) == 700.0d0))
                    passed = (passed .and. (global_sum(2) == 350.0d0))
                    passed = (passed .and. (global_sum(3) == 350.0d0))

                case (4)
                    ! buoyancy per parcel = (1.5 * 2.0 + 0.5 * 1.0)
                    passed = (passed .and. (global_sum(1) == 700.0d0))
                    passed = (passed .and. (global_sum(2) == 350.0d0))
                    passed = (passed .and. (global_sum(3) == 350.0d0))

                case default
                    call mpi_exit_on_error("No valid parcel setup.")
            end select
        end subroutine check_result

end program test_mpi_mixing
