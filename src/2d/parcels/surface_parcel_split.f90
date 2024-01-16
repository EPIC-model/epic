! =============================================================================
!                       Module to split surface parcels
! =============================================================================
module surface_parcel_split
    use options, only : verbose
    use constants, only : f12
    use parameters, only : lmax
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , surface_parcel_reorder         &
                                       , get_surface_parcel_length      &
                                       , surface_parcel_sort
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
            integer                                    :: n, m

            last_index = n_par

            do n = 1, last_index

                ! parcel length:
                l = get_surface_parcel_length(n, sp)

!                 print *, "lmax", lmax, l, (l > lmax), size(sp%vorticity)
!                 stop

                if (l <= lmax) then
                    cycle
                endif

                !
                ! this line is split, i.e., add a new surface parcel
                !

                m = n_par + 1

                ! we only need to add one new parcel
                n_par = m

                sp%volume(n) = f12 * sp%volume(n)

                sp%vorticity(m) = sp%vorticity(n)
                sp%buoyancy(m) = sp%buoyancy(n)
                sp%volume(m) = sp%volume(n)
#ifndef ENABLE_DRY_MODE
                sp%humidity(m) = sp%humidity(n)
#endif
                sp%position(n) = sp%position(n)
                sp%position(m) = sp%position(n) + f12 * l

                sp%right(n) = m
                sp%right(m) = sp%right(n)
            enddo

!             call surface_parcel_reorder(n_par, sp)
            call surface_parcel_sort(n_par, sp)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. surface parcels before and after split: ", last_index, "...", n_par
            endif
#endif
        end subroutine split_lines

end module surface_parcel_split
