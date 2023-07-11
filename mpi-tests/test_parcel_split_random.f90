program test_parcel_split_random
    use mpi_communicator
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    use parcel_container
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_communicate
    use fields, only : field_default
    use parcel_bc, only : apply_periodic_bc
    use parcel_interpl, only : par2grid
    use parcel_split_mod, only : parcel_split
    use mpi_layout, only : box, mpi_layout_init
    use test_utils
    implicit none

    integer, parameter   :: nt = 100
    integer              :: i, n, sk, j, n_orig, n_splits
    integer, allocatable :: seed(:)
    double precision     :: rn(3)
    double precision     :: vol, b(5)

    !--------------------------------------------------------------------------
    ! Initialise MPI and setup all timers:

    call mpi_comm_initialise

    if (comm%rank == comm%master) then
        print '(a35, i6, a11)', "Running 'test_parcel_split_random' with ", comm%size, " MPI ranks."
    endif

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = comm%rank
    call random_seed(put=seed)


    call register_all_timers

    call start_timer(epic_timer)

    !--------------------------------------------------------------------------
    ! Set model parameters:

    nx = 128
    ny = 128
    nz = 64
    lower = (/zero, zero, zero/)
    extent = (/two, two, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    parcel%n_per_cell = 8
    parcel%lambda_max = 4.0d0

    !--------------------------------------------------------------------------
    ! Setup fields: All fields are zero

    call field_default

    !--------------------------------------------------------------------------
    ! Setup parcels:

    call parcel_default

    vol = parcels%volume(1)
    b = parcels%B(:, 1)

    !--------------------------------------------------------------------------
    ! Check initial values:
    call perform_checks

    !--------------------------------------------------------------------------
    ! Do time loop:
    do i = 1, nt

        if (comm%rank == comm%master) then
            print '(a15, i4)', "Performing step", i
        endif

        ! Move each parcel by random value in x and y
        do n = 1, n_parcels
            call random_number(rn)

            if (rn(3) > 0.5d0) then
                call random_number(rn(3))
                j = nint(n_parcels * rn(3)) + 1
                parcels%B(1, j) = 5.0d0 * parcels%B(1, j)
                parcels%buoyancy(j) = 1.0d0
            endif

        enddo

        n_splits = int(sum(parcels%buoyancy(1:n_parcels)))
        call perform_integer_reduction(n_splits)

        if (comm%rank == comm%master) then
            print *, "Split", n_splits, "of", n_total_parcels, "parcels."
        endif

        n_orig = n_parcels

        ! Split parcels
        call parcel_split

        ! Interpolate parcel data to grid
        call par2grid

        ! Test number of parcels: etc.
        call perform_checks

        ! Reset:
        n_parcels = n_orig
        do n = 1, n_parcels
            parcels%volume(n) = vol
            parcels%buoyancy(n) = 0.0d0
            parcels%B(:, n) = b

            call random_number(rn)
            parcels%position(:, n) = box%lower + box%extent * rn
        enddo

        call perform_integer_reduction(n_orig)
        n_total_parcels = n_orig

        ! Do halo swap
        call parcel_communicate

        ! Do periodic shift in x and y
        do n = 1, n_parcels
            call apply_periodic_bc(parcels%position(:, n))
        enddo
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

end program test_parcel_split_random
