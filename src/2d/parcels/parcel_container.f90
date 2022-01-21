! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, hli, center
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
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_delx(x1, x2) result (delx)
#else
        elemental function get_delx(x1, x2) result (delx)
#endif
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
#ifndef NDEBUG
            if ((x1 < lower(1)) .or. (x2 < lower(1)) .or. (x1 > upper(1)) .or. (x2 > upper(1))) then
                write(*,*) 'point outside domain was fed into get_delx'
                write(*,*) 'x1, x2, lower(1), upper(1)'
                write(*,*) x1, x2, lower(1), upper(1)
            endif
#endif
            ! works across periodic edge
            delx = delx - extent(1) * dble(int(delx * hli(1)))
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

            parcels%position(:, n) = parcels%position(:, m)

            parcels%vorticity(n) = parcels%vorticity(m)

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n) = parcels%humidity(m)
#endif
            parcels%B(:, n) = parcels%B(: m)

        end subroutine parcel_replace

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(2, num))
            allocate(parcels%vorticity(num))
            allocate(parcels%B(2, num))
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
