! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, hli, center
    use parcel_ellipse, only : get_angles
    implicit none

    integer :: n_parcels

    type parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            velocity,   &
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12

        double precision, allocatable, dimension(:) :: &
            volume,     &
            vorticity,  &
            buoyancy,   &
            humidity,   &
            stretch
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
            ! works across periodic edge
            delx = delx - extent(1) * dble(int((delx - center(1)) * hli(1)))
        end function get_delx


        ! overwrite parcel n with parcel m
        subroutine parcel_replace(n, m)
            integer, intent(in) :: n, m

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print '(a19, i0, a6, i0)', '    replace parcel ', n, ' with ', m
            endif
#endif

            parcels%position(n, 1) = parcels%position(m, 1)
            parcels%position(n, 2) = parcels%position(m, 2)

            parcels%velocity(n, 1) = parcels%velocity(m, 1)
            parcels%velocity(n, 2) = parcels%velocity(m, 2)

            parcels%vorticity(n) = parcels%vorticity(m)

            if (allocated(parcels%stretch)) then
                parcels%stretch(n)  = parcels%stretch(m)
            endif

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
            parcels%humidity(n) = parcels%humidity(m)

            if (allocated(parcels%B)) then
                parcels%B(n, 1) = parcels%B(m, 1)
                parcels%B(n, 2) = parcels%B(m, 2)
            endif

        end subroutine parcel_replace


        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(num, 2))
            allocate(parcels%velocity(num, 2))
            allocate(parcels%vorticity(num))
            allocate(parcels%stretch(num))
            allocate(parcels%B(num, 2))
            allocate(parcels%volume(num))
            allocate(parcels%buoyancy(num))
            allocate(parcels%humidity(num))
        end subroutine parcel_alloc


        subroutine parcel_dealloc
            deallocate(parcels%position)
            deallocate(parcels%velocity)
            deallocate(parcels%vorticity)

            if (allocated(parcels%stretch)) then
                deallocate(parcels%stretch)
            endif

            if (allocated(parcels%B)) then
                deallocate(parcels%B)
            endif

            deallocate(parcels%volume)
            deallocate(parcels%buoyancy)
            deallocate(parcels%humidity)
        end subroutine parcel_dealloc

end module parcel_container
