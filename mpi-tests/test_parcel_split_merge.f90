program test_parcel_spli_merge
    use mpi_environment
    use options, only : parcel
    use constants, only : zero, one, two
    use parameters, only : update_parameters, nx, ny, nz, lower, extent, vmin, dx
    use parcel_container
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_communicate
    use fields, only : field_default
    use parcel_bc, only : apply_periodic_bc, apply_reflective_bc
    use parcel_interpl, only : par2grid
    use parcel_split_mod, only : parcel_split
    use parcel_merging, only : parcel_merge
    use parcel_nearest
    use mpi_layout, only : mpi_layout_init
    use test_utils
    implicit none

    integer, parameter   :: nt = 100
    integer              :: i, n, sk, n_orig, n_merges
    integer, allocatable :: seed(:)
    double precision     :: rn(3)

    !--------------------------------------------------------------------------
    ! Initialise MPI and setup all timers:

    call mpi_env_initialise

    if (world%rank == world%root) then
        print '(a35, i6, a11)', "Running 'test_parcel_spli_merge' with ", world%size, " MPI ranks."
    endif

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = world%rank
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

    call nearest_win_allocate

    !--------------------------------------------------------------------------
    ! Setup fields: All fields are zero

    call field_default

    !--------------------------------------------------------------------------
    ! Setup parcels:

    call parcel_default

    !--------------------------------------------------------------------------
    ! Check initial values:
    call perform_checks

    !--------------------------------------------------------------------------
    ! Do time loop:
    do i = 1, nt

        if (world%rank == world%root) then
            print '(a15, i4)', "Performing step", i
        endif

        ! Move each parcel by random value in x and y
        do n = 1, n_parcels
            call random_number(rn)
            parcels%position(1:2, n) = parcels%position(1:2, n) + rn(1:2) * dx(1:2)
        enddo

        call parcel_communicate

        do n = 1, n_parcels
            call apply_periodic_bc(parcels%position(:, n))
            call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
        enddo

        n_merges = count(parcels%volume(1:n_parcels) < vmin)
        call perform_integer_reduction(n_merges)

        if (world%rank == world%root) then
            print *, "Merge", n_merges, "of", n_total_parcels, "parcels."
        endif

        call parcel_merge

        do n = 1, n_parcels
            call random_number(rn)
            parcels%B(1, n) = (1.0d0 + rn(1) * 0.3d0) * parcels%B(1, n)
            parcels%B(3, n) = (1.0d0 - rn(2) * 0.3d0) * parcels%B(3, n)
        enddo

        ! Split parcels
        n_orig = n_total_parcels
        call parcel_split

        if (world%rank == world%root) then
            print *, "Split", n_total_parcels - n_orig, "of", n_total_parcels, "parcels."
        endif

        ! Interpolate parcel data to grid
        call par2grid

        ! Test number of parcels: etc.
        call perform_checks
    enddo

    !--------------------------------------------------------------------------
    ! Finish: Free memory and finalise MPI

    call parcel_dealloc

    call nearest_win_deallocate

    call stop_timer(epic_timer)

    call print_timer

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

end program test_parcel_spli_merge
