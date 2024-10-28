! =============================================================================
!                       Test nearest algorithm
!
!               This unit test checks (a, b) - c - d = e and
!               across MPI boundaries for a subcommunicator.
!
! To plot the initial condition use:
! import numpy as np
! import matplotlib.pyplot as plt
!
! x, y, n, r = np.loadtxt('initial.dat', unpack=True)
!
! ng = 16
!
! h = 1.0 / ng
!
! for i in range(ng+1):
!     plt.axhline(i*h, color='k', linestyle='dashed', linewidth=0.5)
!     plt.axvline(i*h, color='k', linestyle='dashed', linewidth=0.5)
!
! colors = ['red', 'blue', 'orange', 'green', 'cyan']
!
! for i, txt in enumerate(n):
!     plt.scatter(x[i], y[i], s=3, color=colors[int(r[i])])
!
! plt.savefig('initial.png', dpi=200)
! plt.close()
!
!
! To plot the final result use:
! import numpy as np
! import matplotlib.pyplot as plt
!
! rank, is_x, is_y, iss, ic_x, ic_y, icc = np.loadtxt('final.dat', unpack=True)
!
! ng = 16
!
! h = 1.0 / ng
! colors = ['red', 'blue', 'orange', 'green', 'cyan', 'magenta']
! marker = ['o', '+', '*', '.', '.', '.']
!
!
! for i in range(ng+1):
!     plt.axhline(i*h, color='k', linestyle='dashed', linewidth=0.5)
!     plt.axvline(i*h, color='k', linestyle='dashed', linewidth=0.5)
!
! for i, r in enumerate(rank):
!     plt.scatter(is_x[i], is_y[i], marker='o', color=colors[int(r)], s=15, alpha=0.25)
!     plt.scatter(ic_x[i], ic_y[i], marker='x', color=colors[int(r)], s=15)
! plt.savefig('final.png', dpi=500)
! plt.close()
! =============================================================================
program test_mpi_nearest_subcomm
    use unit_test
    use constants, only : pi, zero, one, two, five, ten
    use parcels_mod, only : parcels
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, dx, vmin, max_num_parcels
    use parcel_nearest
    use mpi_environment
    use mpi_layout
    use mpi_timer
    use mpi_utils, only : mpi_exit_on_error
    use mpi_utils, only : mpi_stop
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma, inva
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, check_array(2), n_invalid
!     integer :: is, ic, r

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 16
    ny = 16
    nz = 16
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

    !
    ! Test 1: One patch in the centre (a, b) - c - d = e
    !
    call parcel_setup(2)

    call find_nearest(parcels, isma, iclo, inva, n_merge, n_invalid)

!     ! To print out the result enable the following lines:
!     do r = 0, world%size-1
!         if (r == world%rank) then
!             do n = 1, n_merge
!                 is = isma(n)
!                 ic = iclo(n)
!                 print *, world%rank, parcels%position(1, is), parcels%position(2, is), int(parcels%buoyancy(is)), &
!                                     parcels%position(1, ic), parcels%position(2, ic), int(parcels%buoyancy(ic))
!             enddo
!         endif
!         call MPI_Barrier(world%comm, world%err)
!     enddo
!     call mpi_stop

    call check_result(300)

    if (world%rank == world%root) then
        call print_result_logical('Test MPI nearest algorithm subcomm', passed)
    endif

    call nearest_win_deallocate

    call mpi_env_finalise

    contains

        ! One patch (4x4 grid cells) of small parcels in the centre:
        ! (a, b) - c - d = e
        subroutine cell_placement_1(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z, fac

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            fac = 1.2d0

            if ((i > 5) .and. (i < 10) .and. (j > 5) .and. (j < 10)) then
                fac = 0.9d0
            endif

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.4d0
            parcels%position(2, l) = y - dx(2) * 0.35d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x - dx(1) * 0.35d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel e
            parcels%position(1, l) = x - dx(1) * 0.27d0
            parcels%position(2, l) = y + dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_1

        ! Multiple patches (4x4 grid cells) of small parcels:
        ! (a, b) - c - d = e
        subroutine cell_placement_2(l, i, j, k)
            integer, intent(inout) :: l
            integer, intent(in)    :: i, j, k
            integer                :: ix, iy, iz
            double precision       :: x, y, z, fac

            ix = i
            iy = j
            iz = k

            x = lower(1) + (0.5d0 + dble(ix)) * dx(1)
            y = lower(2) + (0.5d0 + dble(iy)) * dx(2)
            z = lower(3) + (0.5d0 + dble(iz)) * dx(3)

            fac = 1.2d0

            if ((i >= 0) .and. (i < 4) .and. (j >= 0) .and. (j < 4)) then
                fac = 0.9d0
            else if ((i >= 9) .and. (i < 14) .and. (j >= 3) .and. (j < 7)) then
                fac = 0.9d0
            else if ((i > 2) .and. (i < 6) .and. (j > 10) .and. (j < 14)) then
                fac = 0.9d0
            else
                return
            endif

            ! small parcel a
            parcels%position(1, l) = x + dx(1) * 0.4d0
            parcels%position(2, l) = y - dx(2) * 0.35d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel b
            parcels%position(1, l) = x + dx(1) * 0.44d0
            parcels%position(2, l) = y + dx(2) * 0.4d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel c
            parcels%position(1, l) = x - dx(1) * 0.45d0
            parcels%position(2, l) = y - dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel d
            parcels%position(1, l) = x - dx(1) * 0.35d0
            parcels%position(2, l) = y + dx(2) * 0.44d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

            ! small parcel e
            parcels%position(1, l) = x - dx(1) * 0.27d0
            parcels%position(2, l) = y + dx(2) * 0.38d0
            parcels%position(3, l) = z
            parcels%volume(l) = fac * vmin
            parcels%buoyancy(l) = l + world%rank * 100
            l = l + 1

        end subroutine cell_placement_2

        subroutine parcel_setup(num)
            integer, intent(in) :: num
            integer             :: i, j, k
!             integer             :: r

            n = 1
            do k = box%lo(3)+1, box%lo(3)+1
                do j = box%lo(2), box%hi(2)
                    do i = box%lo(1), box%hi(1)
                        select case (num)
                            case (1)
                                call cell_placement_1(n, i, j, k)
                            case (2)
                                call cell_placement_2(n, i, j, k)
                            case default
                                call mpi_exit_on_error("No valid parcel setup.")
                        end select
                    enddo
                enddo
            enddo

!             ! To print out the parcel setups enable the following lines:
!             do r = 0, world%size-1
!                 if (r == world%rank) then
!                     do i = 1, n - 1
!                         print *, parcels%position(1, i), parcels%position(2, i), i + 100 * world%rank, world%rank
!                     enddo
!                 endif
!                 call MPI_Barrier(world%comm, world%err)
!             enddo
!             call MPI_Barrier(world%comm, world%err)
!             call mpi_stop

            parcels%local_num = n - 1
            parcels%total_num = 0

            call MPI_Allreduce(parcels%local_num,   &
                               parcels%total_num,   &
                               1,                   &
                               MPI_INTEGER,         &
                               MPI_SUM,             &
                               world%comm,          &
                               world%err)
        end subroutine parcel_setup

        subroutine check_result(n_true_merges)
            integer, intent(in) :: n_true_merges
            check_array(1) = parcels%local_num - n_invalid
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
                passed = (passed .and. (check_array(1) == parcels%total_num) .and. (check_array(2) == n_true_merges))
            endif
        end subroutine check_result

end program test_mpi_nearest_subcomm
