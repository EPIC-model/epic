! =============================================================================
!                       Module to split surface parcels
! =============================================================================
module surface_parcel_split
    use options, only : verbose
    use constants, only : pi, three, f12, f14
    use parameters, only : lmax
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , surface_parcel_reorder
    use omp_lib
    implicit none

    contains

        ! Split large parcels (volumes larger than lmax) or
        ! parcels with aspect ratios larger than the threshold.
        ! @param[inout] parcels
        subroutine split_lines(n_par, sp)
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: sp
            double precision                           :: l
            integer                                    :: last_index
            integer                                    :: n, m, j

            last_index = n_par

            do n = 1, last_index

                j = sp%right(n)

                ! parcel length:
                l = sp%position(j) - sp%position(n)

                if (l <= lmax) then
                    cycle
                endif

                !
                ! this lines is split, i.e., add a new surface parcel
                !

                m = n_par + 1

                ! we only need to add one new parcel
                n_par = n_par + 1

                sp%vorticity(m) = sp%vorticity(n)
                sp%buoyancy(m) = sp%buoyancy(n)
#ifndef ENABLE_DRY_MODE
                sp%humidity(m) = sp%humidity(n)
#endif
                sp%position(n) = sp%position(n)
                sp%position(m) = sp%position(n) + f12 * l

                sp%right(n) = m
                sp%right(m) = j
            enddo

            call surface_parcel_reorder(n_par, sp)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. surface parcels before and after split: ", last_index, "...", n_par
            endif
#endif
        end subroutine split_lines

end module surface_parcel_split
