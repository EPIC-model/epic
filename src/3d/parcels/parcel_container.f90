! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, extenti, hl
    implicit none

    integer :: n_parcels

    type parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            vorticity,  &
            B               ! B matrix entries; ordering:
                            ! B(:, 1) = B11, B(:, 2) = B12, B(:, 3) = B13
                            ! B(:, 4) = B22, B(:, 5) = B23

        double precision, allocatable, dimension(:) :: &
            volume,     &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        ! Obtain the difference between two zonal coordinates
        ! across periodic edges
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
        ! @returns delx = x1 - x2
        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
            ! works across periodic edge
            delx = delx - extent(1) * dble(floor((delx+hl(1))*extenti(1)))
        end function get_delx

        ! Obtain the difference between two meridional coordinates
        ! across periodic edges
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        elemental function get_dely(y1, y2) result (dely)
            double precision, intent(in) :: y1, y2
            double precision             :: dely

            dely = y1 - y2
            ! works across periodic edge
            dely = dely - extent(2) * dble(floor((dely+hl(2))*extenti(2)))
        end function get_dely


        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine parcel_replace(n, m)
            integer, intent(in) :: n, m

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print '(a19, i0, a6, i0)', '    replace parcel ', n, ' with ', m
            endif
#endif

            parcels%position(:, n) = parcels%position(:, m)

            parcels%vorticity(:, n) = parcels%vorticity(:, m)

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n) = parcels%humidity(m)
#endif
            parcels%B(:, n) = parcels%B(:, m)

        end subroutine parcel_replace

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(3, num))
            allocate(parcels%vorticity(3, num))
            allocate(parcels%B(5, num))
            allocate(parcels%volume(num))
            allocate(parcels%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(parcels%humidity(num))
#endif
        end subroutine parcel_alloc

        ! Deallocate parcel memory
        subroutine parcel_dealloc
            deallocate(parcels%position)
            deallocate(parcels%vorticity)
            deallocate(parcels%B)
            deallocate(parcels%volume)
            deallocate(parcels%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(parcels%humidity)
#endif
        end subroutine parcel_dealloc

end module parcel_container
