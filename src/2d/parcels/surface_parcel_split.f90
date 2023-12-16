! =============================================================================
!                       Module to split surface parcels
! =============================================================================
module surface_parcel_split
    use options, only : verbose
    use constants, only : pi, three, f12, f14
    use parameters, only : lmax
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , surface_parcel_sort
    use omp_lib
    implicit none

    contains

        ! Split large parcels (volumes larger than lmax) or
        ! parcels with aspect ratios larger than the threshold.
        ! @param[inout] parcels
        subroutine split_lines(n_par, parcels)
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: parcels
            double precision                           :: l
            integer                                    :: last_index
            integer                                    :: n, m, j

            last_index = n_par

            do n = 1, last_index

                j = parcels%right(n)

                ! parcel length:
                l = parcels%position(j) - parcels%position(n)

                if (l <= lmax) then
                    cycle
                endif

                !
                ! this lines is split, i.e., add a new surface parcel
                !

                m = n_par + 1

                ! we only need to add one new parcel
                n_par = n_par + 1

                parcels%vorticity(m) = parcels%vorticity(n)
                parcels%buoyancy(m) = parcels%buoyancy(n)
#ifndef ENABLE_DRY_MODE
                parcels%humidity(m) = parcels%humidity(n)
#endif
                parcels%position(n) = parcels%position(n)
                parcels%position(m) = parcels%position(n) + f12 * l

                parcels%right(n) = m
                parcels%right(m) = j
            enddo

            call surface_parcel_sort(n_par, parcels)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. surface parcels before and after split: ", last_index, "...", n_par
            endif
#endif
        end subroutine split_lines

end module surface_parcel_split
