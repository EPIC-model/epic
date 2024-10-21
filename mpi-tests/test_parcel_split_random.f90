program test_parcel_split_random
    use mpi_environment
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent, max_num_parcels
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
    integer              :: i, n, j, n_orig, n_splits
    double precision     :: rn(3)
    double precision     :: vol, b(5)
    logical, allocatable :: picked(:)

    !--------------------------------------------------------------------------
    ! Initialise MPI and setup all timers:

    call mpi_env_initialise

    if (world%rank == world%root) then
        print '(a35, i6, a11)', "Running 'test_parcel_split_random' with ", world%size, " MPI ranks."
    endif

    call init_rng

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

    allocate(picked(max_num_parcels))

    picked = .false.

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

        if (world%rank == world%root) then
            print '(a15, i4)', "Performing step", i
        endif

        do n = 1, n_parcels
            call random_number(rn)

            if (rn(3) > 0.5d0) then
                call random_number(rn(3))
                j = nint(n_parcels * rn(3)) + 1

                if (.not. picked(j)) then
                    parcels%B(1, j) = 5.0d0 * parcels%B(1, j)

                    picked(j) = .true.
                endif
            endif

        enddo

        n_splits = count(picked)
        call perform_integer_reduction(n_splits)

        if (world%rank == world%root) then
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
            parcels%B(:, n) = b

            call random_number(rn)
            parcels%position(:, n) = box%lower + box%extent * rn
        enddo

        picked = .false.

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

    deallocate(picked)

    call mpi_env_finalise


    contains

        subroutine perform_checks
            call check_total_number_of_parcels

        end subroutine perform_checks

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine check_total_number_of_parcels
            integer :: n_total

            n_total = n_parcels

            call perform_integer_reduction(n_total)

            if (world%rank == world%root) then
                if (n_total /= n_total_parcels) then
                    print *, "check_total_number_of_parcels: Total number of parcels disagree!"
                    call MPI_Abort(world%comm, -1, world%err)
                endif
            endif

        end subroutine check_total_number_of_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine perform_integer_reduction(var)
            integer, intent(inout) :: var

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE, var, 1, MPI_INTEGER, MPI_SUM, &
                                world%root, world%comm, world%err)
            else
                call MPI_Reduce(var, var, 1, MPI_INTEGER, MPI_SUM, &
                                world%root, world%comm, world%err)
            endif

        end subroutine perform_integer_reduction

end program test_parcel_split_random
