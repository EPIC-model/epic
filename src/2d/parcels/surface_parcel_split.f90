! =============================================================================
!                       Module to split surface parcels
! =============================================================================
module surface_parcel_split
    use options, only : verbose, parcel
    use constants, only : pi, three, f12, f14
    use parameters, only : lmax
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_top_parcels, top_parcels     &
                                       , n_bot_parcels, bot_parcels
    use omp_lib
    implicit none

    private :: split_lines

    contains

        subroutine split_surface_parcels

            call split_lines(n_top_parcels, top_parcels)
            call split_lines(n_bot_parcels, bot_parcels)

        end subroutine split_surface_parcels

        ! Split large parcels (volumes larger than lmax) or
        ! parcels with aspect ratios larger than the threshold.
        ! @param[inout] parcels
        subroutine split_lines(n_par, parcels)
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: parcels
            double precision                           :: l, h
            integer                                    :: last_index
            integer                                    :: n, n_thread_loc

            last_index = n_par

            !$omp parallel default(shared)
            !$omp do private(n, l, h, n_thread_loc)
            do n = 1, last_index

                l = parcels%length(n)

                if (l <= parcel%lambda_max .and. l <= lmax) then
                    cycle
                endif

                !
                ! this lines is split, i.e., add a new surface parcel
                !

                h = f14 * l
                parcels%length(n) = f12 * l

                !$omp critical
                n_thread_loc = n_par + 1

                ! we only need to add one new parcel
                n_par = n_par + 1
                !$omp end critical

                parcels%vorticity(n_thread_loc) = parcels%vorticity(n)
                parcels%length(n_thread_loc) = parcels%length(n)
                parcels%buoyancy(n_thread_loc) = parcels%buoyancy(n)
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n_thread_loc) = parcels%humidity(n)
#endif
                parcels%position(n_thread_loc) = parcels%position(n) - h
                parcels%position(n) = parcels%position(n) + h
            enddo
            !$omp end do
            !$omp end parallel

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. surface parcels before and after split: ", last_index, "...", n_par
            endif
#endif
        end subroutine split_lines

end module surface_parcel_split
