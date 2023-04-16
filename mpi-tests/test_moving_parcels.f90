program test_moving_parcels
    use mpi_communicator
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent, dx
    use parcel_container
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_halo_swap
    use fields, only : field_default
    use parcel_bc, only : apply_periodic_bc
    use test_utils
    implicit none

    integer, parameter :: nt = 100
    integer            :: i, n

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


    !--------------------------------------------------------------------------
    ! Do time loop:

    do i = 1, nt

        if (comm%rank == comm%master) then
            print '(a15, i4)', "Performing step", i
        endif

        ! Move each parcel by dx and dy
        do n = 1, n_parcels
            parcels%position(1:2, n) = parcels%position(1:2, n) + dx(1:2)
        enddo

        ! Do halo swap
        call parcel_halo_swap

        ! Do periodic shift in x and y
        do n = 1, n_parcels
            call apply_periodic_bc(parcels%position(:, n))
        enddo

        ! Test number of parcels: etc.
        call check_total_number_of_parcels
    enddo

    !--------------------------------------------------------------------------
    ! Finish: Free memory and finalise MPI

    call parcel_dealloc

    call stop_timer(epic_timer)

    call print_timer

    call mpi_comm_finalise


    contains

        subroutine check_total_number_of_parcels
            integer :: n_total

            n_total = n_parcels

            if (comm%rank == comm%master) then
                call MPI_Reduce(MPI_IN_PLACE, n_total, 1, MPI_INTEGER, MPI_SUM, &
                                comm%master, comm%world, comm%err)
            else
                call MPI_Reduce(n_total, n_parcels, 1, MPI_INTEGER, MPI_SUM, &
                                comm%master, comm%world, comm%err)
            endif

            if (comm%rank == comm%master) then
                if (n_total /= n_total_parcels) then
                    print *, "Total number of parcels disagree!"
                    call MPI_Abort(comm%world, -1, comm%err)
                endif
            endif

        end subroutine check_total_number_of_parcels

end program test_moving_parcels
