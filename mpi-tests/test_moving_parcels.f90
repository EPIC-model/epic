program test_moving_parcels
    use mpi_communicator
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent, dx
    use parcel_container
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_halo_swap
    use fields, only : field_default, nparg, vortg
    use parcel_bc, only : apply_periodic_bc
    use mpi_layout, only : box
    use parcel_interpl, only : par2grid
    use test_utils
    implicit none

    integer, parameter :: nt = 100
    integer            :: i, n
    double precision   :: VOR(3) = (/1.5d0, 2.0d0, 2.5d0/)

    !--------------------------------------------------------------------------
    ! Initialise MPI and setup all timers:

    call mpi_comm_initialise

    if (comm%rank == comm%master) then
        print '(a35, i6, a11)', "Running test 'moving parcels' with ", comm%size, " MPI ranks."
    endif

    call register_all_timers

    call start_timer(epic_timer)

    !--------------------------------------------------------------------------
    ! Set model parameters:

    nx = 128
    ny = 128
    nz = 64
    lower = (/zero, zero, zero/)
    extent = (/two, two, one/)

    call update_parameters

    parcel%n_per_cell = 8

    !--------------------------------------------------------------------------
    ! Setup fields: All fields are zero

    call field_default

    !--------------------------------------------------------------------------
    ! Setup parcels:

    call parcel_default

    do i = 1, 3
        parcels%vorticity(i, 1:n_parcels) = VOR(i)
    enddo

    !--------------------------------------------------------------------------
    ! Interpolate parcel data to grid:
    call par2grid

    !--------------------------------------------------------------------------
    ! Check initial values:
    call perform_checks

    !--------------------------------------------------------------------------
    ! Do time loop:

    do i = 1, nt

        if (comm%rank == comm%master) then
            print '(a15, i4)', "Performing step", i
        endif

        VOR = 1.00d0 * VOR

        ! Move each parcel by dx and dy
        do n = 1, n_parcels
            parcels%position(1:2, n) = parcels%position(1:2, n) + dx(1:2)
            parcels%vorticity(:, n) = VOR
        enddo

        ! Do halo swap
        call parcel_halo_swap

        ! Do periodic shift in x and y
        do n = 1, n_parcels
            call apply_periodic_bc(parcels%position(:, n))
        enddo

        ! Interpolate parcel data to grid
        call par2grid

        ! Test number of parcels: etc.
        call perform_checks
    enddo

    !--------------------------------------------------------------------------
    ! Finish: Free memory and finalise MPI

    call parcel_dealloc

    call stop_timer(epic_timer)

    call print_timer

    call mpi_comm_finalise


    contains

        subroutine perform_checks
            call check_total_number_of_parcels

            call check_number_of_parcels_per_cell

            call check_gridded_vorticity

        end subroutine perform_checks

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine check_total_number_of_parcels
            integer :: n_total

            n_total = n_parcels

            call perform_integer_reduction(n_total)

            if (comm%rank == comm%master) then
                if (n_total /= n_total_parcels) then
                    print *, "check_total_number_of_parcels: Total number of parcels disagree!"
                    call MPI_Abort(comm%world, -1, comm%err)
                endif
            endif

        end subroutine check_total_number_of_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine check_number_of_parcels_per_cell
            integer :: n_total

            n_total = sum(nparg(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            call perform_integer_reduction(n_total)

            if (comm%rank == comm%master) then
                if (n_total /= n_total_parcels) then
                    print *, "check_number_of_parcels_per_cell: Total number of parcel disagrees with grid!"
                    call MPI_Abort(comm%world, -1, comm%err)
                endif
            endif

        end subroutine check_number_of_parcels_per_cell

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine check_gridded_vorticity
            double precision :: vmin, vmax
            integer          :: j
            character(1)     :: dir(3) = (/'x', 'y', 'z'/)

            do j = 1, 3
                vmin = minval(vortg(0:nz, :, :, j))
                vmax = maxval(vortg(0:nz, :, :, j))

                if ((dabs(vmin - VOR(j)) > 100.0d0 * epsilon(vmin)) .or. &
                    (dabs(vmax - VOR(j)) > 100.0d0 * epsilon(vmax))) then
                    print *, "check_gridded_vorticity: Error in " // dir(j) // "-vorticity!"
                    call MPI_Abort(comm%world, -1, comm%err)
                endif
            enddo
        end subroutine check_gridded_vorticity

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine perform_integer_reduction(var)
            integer, intent(inout) :: var

            if (comm%rank == comm%master) then
                call MPI_Reduce(MPI_IN_PLACE, var, 1, MPI_INTEGER, MPI_SUM, &
                                comm%master, comm%world, comm%err)
            else
                call MPI_Reduce(var, var, 1, MPI_INTEGER, MPI_SUM, &
                                comm%master, comm%world, comm%err)
            endif

        end subroutine perform_integer_reduction

end program test_moving_parcels
