! =============================================================================
!                        Submodule to split ellipses
! =============================================================================
submodule (parcel_ellipse) parcel_ellipse_split_smod
    use options, only : verbose, parcel
    use constants, only : three, f34
    use parameters, only : amax
    use mpi_environment
    use mpi_collectives
!     use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

!     integer :: surf_split_timer

    contains

        ! Split large parcels (areas larger than amax) or
        ! parcels with aspect ratios larger than the threshold.
        ! @param[inout] parcels
        ! @param[in] threshold is the largest allowed aspect ratio
        subroutine parcel_ellipse_split(this)
            class(ellipse_pc_type), intent(inout) :: this
            double precision                      :: a2, lam
            double precision                      :: evec(2)
            double precision                      :: h
            integer                               :: last_index
            integer                               :: n, n_thread_loc
            integer                               :: pid(2 * this%local_num)
            integer, allocatable                  :: invalid(:)

!             call start_timer(surf_split_timer)


            last_index = this%local_num

            !$omp parallel default(shared)
            !$omp do private(n, a2, lam, evec, h, n_thread_loc)
            do n = 1, last_index

                a2 = this%get_eigenvalue(n)

                ! a/b
                lam = this%get_aspect_ratio(n, a2)

                if (lam <= parcel%lambda_max .and. this%area(n) <= amax) then
                    cycle
                endif

                pid(n) = 0

                !
                ! this ellipse is split, i.e., add a new parcel
                !

                evec = this%get_eigenvector(n, a2)

                this%B(1, n) = this%B(1, n) - f34 * a2 * evec(1) ** 2
                this%B(2, n) = this%B(2, n) - f34 * a2 * (evec(1) * evec(2))
                this%B(3, n) = this%B(3, n) - f34 * a2 * evec(2) ** 2

                h = f14 * dsqrt(three * a2)
                this%area(n) = f12 * this%area(n)

                !$omp critical
                n_thread_loc = this%local_num + 1

                ! we only need to add one new parcel
                this%local_num = this%local_num + 1
                !$omp end critical


                this%B(:, n_thread_loc) = this%B(:, n)

                this%area(n_thread_loc) = this%area(n)
                this%position(:, n_thread_loc) = this%position(:, n) - h * evec
                this%position(:, n) = this%position(:, n) + h * evec

                ! save parcel indices of child parcels for the
                ! halo swap routine
                pid(n) = n
                pid(n_thread_loc) = n_thread_loc
            enddo
            !$omp end do
            !$omp end parallel

            ! after this operation the root MPI process knows the new
            ! number of parcels in the simulation
            this%total_num = this%local_num
            call mpi_blocking_reduce(this%total_num, MPI_SUM, world)

            ! all entries in "pid" that are non-zero are indices of
            ! child parcels; remove all zero entries such that
            ! we can do a halo swap
            invalid = pack(pid(1:this%local_num), pid(1:this%local_num) /= 0)

            ! send the invalid parcels to the proper MPI process;
            ! delete them on *this* MPI process and
            ! apply periodic boundary condition
!             call parcel_communicate(invalid)  !FIXME

!             call stop_timer(surf_split_timer)

        end subroutine parcel_ellipse_split

end submodule parcel_ellipse_split_smod
