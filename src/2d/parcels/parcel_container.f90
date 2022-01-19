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
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12

        double precision, allocatable, dimension(:) :: &
            volume,     &
            vorticity,  &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        ! Obtain the difference between two horizontal coordinates
        ! across periodic edges
        ! @param[in] x1 first horizontal position
        ! @param[in] x2 second horizontal position
        ! @returns delx = x1 - x2
        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
            ! works across periodic edge
            delx = delx - extent(1) * dble(floor((delx+hl(1)) * extenti(1)))

        end function get_delx


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

            parcels%position(n, 1) = parcels%position(m, 1)
            parcels%position(n, 2) = parcels%position(m, 2)

            parcels%vorticity(n) = parcels%vorticity(m)

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n) = parcels%humidity(m)
#endif
            parcels%B(n, 1) = parcels%B(m, 1)
            parcels%B(n, 2) = parcels%B(m, 2)

        end subroutine parcel_replace

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(num, 2))
            allocate(parcels%vorticity(num))
            allocate(parcels%B(num, 2))
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
