! =============================================================================
!                           Module to split ellipses
! =============================================================================
module surface_parcel_split_mod
    use options, only : verbose, parcel
    use constants, only : pi, three, f12, f14, f34
    use parameters, only : amax
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_lo_surf_parcels              &
                                       , n_up_surf_parcels              &
                                       , lo_surf_parcels                &
                                       , up_surf_parcels
    use parcel_ellipse, only : get_eigenvalue      &
                             , get_eigenvector     &
                             , get_aspect_ratio
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: surf_split_timer

    contains

        subroutine surface_parcel_split
            call do_ellipse_split(lo_surf_parcels, n_lo_surf_parcels)
            call do_ellipse_split(up_surf_parcels, n_up_surf_parcels)
        end subroutine surface_parcel_split

        ! Split large parcels (areas larger than amax) or
        ! parcels with aspect ratios larger than parcel%lambda_max.
        ! @param[inout] parcels
        subroutine do_ellipse_split(s_parcels, n_par)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            double precision                                   :: B(3)
            double precision                                   :: a2, lam, V
            double precision                                   :: evec(2)
            double precision                                   :: h
            integer                                            :: last_index
            integer                                            :: n, n_thread_loc

            call start_timer(surf_split_timer)

            last_index = n_par

            !$omp parallel default(shared)
            !$omp do private(n, B, a2, lam, V, evec, h, n_thread_loc)
            do n = 1, last_index
                B(:) = s_parcels%B(:, n)
                V = s_parcels%area(n)

                a2 = get_eigenvalue(B)

                ! a/b
                lam = get_aspect_ratio(a2, V)

                if (lam <= parcel%lambda_max .and. V <= amax) then
                    cycle
                endif

                !
                ! this ellipse is split, i.e., add a new parcel
                !

                evec = get_eigenvector(a2, B)

                s_parcels%B(1, n) = B(1) - f34 * a2 * evec(1) ** 2
                s_parcels%B(2, n) = B(2) - f34 * a2 * (evec(1) * evec(2))
                s_parcels%B(3, n) = B(3) - f34 * a2 * evec(2) ** 2

                h = f14 * dsqrt(three * a2)
                s_parcels%area(n) = f12 * V

                !$omp critical
                n_thread_loc = n_par + 1

                ! we only need to add one new parcel
                n_par = n_par + 1
                !$omp end critical


                s_parcels%B(:, n_thread_loc) = s_parcels%B(:, n)

                s_parcels%area(n_thread_loc) = s_parcels%area(n)
                s_parcels%position(:, n_thread_loc) = s_parcels%position(:, n) - h * evec
                s_parcels%position(:, n) = s_parcels%position(:, n) + h * evec
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(surf_split_timer)

        end subroutine do_ellipse_split

end module surface_parcel_split_mod
