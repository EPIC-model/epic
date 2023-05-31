! =============================================================================
!                           Module to split ellipsoids
! =============================================================================
module parcel_split_mod
    use options, only : verbose, parcel
    use constants, only : pi, three, five, f12, f34
    use parameters, only : vmax, max_num_parcels
    use parcel_container, only : parcels                &
                               , n_parcels              &
                               , n_total_parcels        &
                               , parcel_resize
    use parcel_bc, only : apply_reflective_bc, apply_swap_periodicity
    use parcel_ellipsoid, only : diagonalise, get_aspect_ratio
    use mpi_timer, only : start_timer, stop_timer
    use omp_lib
    use mpi_communicator, only : comm, MPI_SUM
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    double precision, parameter :: dh = f12 * dsqrt(three / five)

    ! saves the last previous number of splits
    ! to calculate the "running average"
    ! (see https://en.wikipedia.org/wiki/Moving_average)
    integer, protected :: n_previous_splits(10) = 0

    ! current index for "running average"
    integer, protected :: ri = 0

    private :: dh, n_previous_splits, ri

    integer :: split_timer

    ! number of parcel splits (is reset in every write step)
    integer :: n_parcel_splits = 0


    contains

        ! Split large parcels (volumes larger than vmax) or
        ! parcels with aspect ratios larger than parcel%lambda_max.
        subroutine parcel_split
            double precision     :: B(5)
            double precision     :: vol, lam
            double precision     :: D(3), V(3, 3)
            integer              :: last_index
            integer              :: n, n_thread_loc
            integer              :: pid(2 * n_parcels)
            integer, allocatable :: invalid(:)
            integer              :: n_estimate, shrunk_size, grown_size
#ifdef ENABLE_VERBOSE
            integer              :: orig_num

            orig_num = n_total_parcels
#endif

            call start_timer(split_timer)

            ! Estimate final number of parcels based on previous number of splits:
            ! If the est
            n_estimate = nint(dble(sum(n_previous_splits)) / dble(size(n_previous_splits))) &
                       + n_parcels

            shrunk_size = nint(parcel%shrink_factor * max_num_parcels)

            if (n_estimate > max_num_parcels) then
                grown_size = nint(parcel%grow_factor * max_num_parcels)
                call parcel_resize(grown_size)
            else if (n_estimate < shrunk_size) then
                call parcel_resize(shrunk_size)
            endif

            last_index = n_parcels

            !$omp parallel default(shared)
            !$omp do private(n, B, vol, lam, D, V, n_thread_loc)
            do n = 1, last_index
                B = parcels%B(:, n)
                vol = parcels%volume(n)

                call diagonalise(B, vol, D, V)

                ! evaluate maximum aspect ratio (a2 >= b2 >= c2)
                lam = get_aspect_ratio(D)

                pid(n) = 0

                if (lam <= parcel%lambda_max .and. vol <= vmax) then
                    cycle
                endif

                !
                ! this ellipsoid is split, i.e., add a new parcel
                !
                parcels%B(1, n) = B(1) - f34 * D(1) * V(1, 1) ** 2
                parcels%B(2, n) = B(2) - f34 * D(1) * V(1, 1) * V(2, 1)
                parcels%B(3, n) = B(3) - f34 * D(1) * V(1, 1) * V(3, 1)
                parcels%B(4, n) = B(4) - f34 * D(1) * V(2, 1) ** 2
                parcels%B(5, n) = B(5) - f34 * D(1) * V(2, 1) * V(3, 1)

                parcels%volume(n) = f12 * vol

                !$omp critical
                n_thread_loc = n_parcels + 1

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1
                !$omp end critical

                parcels%B(:, n_thread_loc) = parcels%B(:, n)

                parcels%vorticity(:, n_thread_loc) = parcels%vorticity(:, n)
                parcels%volume(n_thread_loc) = parcels%volume(n)
                parcels%buoyancy(n_thread_loc) = parcels%buoyancy(n)
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n_thread_loc) = parcels%humidity(n)
#endif
                V(:, 1) = V(:, 1) * dh * dsqrt(D(1))
                parcels%position(:, n_thread_loc) = parcels%position(:, n) - V(:, 1)
                parcels%position(:, n) = parcels%position(:, n) + V(:, 1)

                ! child parcels need to be reflected into domain, if their center
                ! is inside the halo region
                call apply_reflective_bc(parcels%position(:, n_thread_loc), &
                                         parcels%B(:, n_thread_loc))

                call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))

                ! save parcel indices of child parcels for the
                ! halo swap routine
                pid(n) = n
                pid(n_thread_loc) = n_thread_loc
            enddo
            !$omp end do
            !$omp end parallel

            ! number of previous splits; used for calculating
            ! the running average for parcel container size adaption
            n_previous_splits(ri+1) = n_parcels - last_index
            ri = mod(ri + 1, 10)

            n_parcel_splits = n_parcel_splits + n_parcels - last_index

            ! after this operation the root MPI process knows the new
            ! number of parcels in the simulation
            n_total_parcels = n_parcels
            call mpi_blocking_reduce(n_total_parcels, MPI_SUM)

            ! all entries in "pid" that are non-zero are indices of
            ! child parcels; remove all zero entries such that
            ! we can do a halo swap
            invalid = pack(pid(1:n_parcels), pid(1:n_parcels) /= 0)

            ! send the invalid parcels to the proper MPI process;
            ! delete them on *this* MPI process and
            ! apply periodic boundary condition
            call apply_swap_periodicity(invalid)

#ifdef ENABLE_VERBOSE
            if (verbose .and. (comm%rank == comm%master)) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", orig_num, "...", n_total_parcels
            endif
#endif
            call stop_timer(split_timer)

        end subroutine parcel_split

end module parcel_split_mod
