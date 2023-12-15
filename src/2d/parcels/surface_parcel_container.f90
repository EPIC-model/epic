! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module surface_parcel_container
    use options, only : verbose
    use parameters, only : extent, extenti, center, lower, upper
    implicit none

    integer :: n_top_parcels, n_bot_parcels

    type surface_parcel_container_type
        double precision, allocatable, dimension(:) :: &
            position,   &
            length,     &
            vorticity,  &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy
    end type surface_parcel_container_type

    type(surface_parcel_container_type) top_parcels, bot_parcels


    contains

        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine surface_parcel_replace(n, m)
            integer, intent(in) :: n, m

            parcels%position(n) = parcels%position(m)

            parcels%vorticity(n) = parcels%vorticity(m)

            parcels%length(n)  = parcels%length(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n) = parcels%humidity(m)
#endif
        end subroutine surface_parcel_replace

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine surface_parcel_alloc(num)
            integer, intent(in) :: num

            call _alloc(num, top_parcels)
            call _alloc(num, bot_parcels)

        end subroutine surface_parcel_alloc


        ! Deallocate parcel memory
        subroutine surface_parcel_dealloc

            call _dealloc(top_parcels)
            call _dealloc(bot_parcels)

        end subroutine surface_parcel_dealloc

        subroutine _alloc(num, sp)
            integer,                             intent(in)    :: num
            type(surface_parcel_container_type), intent(inout) :: sp

            allocate(sp%position(num))
            allocate(sp%vorticity(num))
            allocate(sp%length(num))
            allocate(sp%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(sp%humidity(num))
#endif
        end subroutine _alloc

        subroutine _dealloc(sp)
            type(surface_parcel_container_type), intent(inout) :: sp

            deallocate(sp%position)
            deallocate(sp%vorticity)
            deallocate(sp%length)
            deallocate(sp%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(sp%humidity)
#endif
        end subroutine _dealloc

end module surface_parcel_container
