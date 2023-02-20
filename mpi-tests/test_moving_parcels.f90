program test_moving_parcels
    use mpi_communicator
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent, dx
    use parcel_container
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_halo_swap
    use fields, only : field_default
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

        ! Test number of parcels: etc.
        ! TODO
    enddo

    !--------------------------------------------------------------------------
    ! Finish: Free memory and finalise MPI

    call parcel_dealloc

    call stop_timer(epic_timer)

    call print_timer

    call mpi_comm_finalise

end program test_moving_parcels
