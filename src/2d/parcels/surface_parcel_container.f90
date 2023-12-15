! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module surface_parcel_container
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

        ! LS-RK-4 arrays
        double precision, allocatable, dimension(:) :: &
            strain,    &   ! strain at parcel location
            delta_pos, &
            delta_vor      ! vorticity integration

    end type surface_parcel_container_type

    type(surface_parcel_container_type) top_parcels, bot_parcels

    private :: alloc, dealloc

    contains

        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine surface_parcel_replace(n, m, sp)
            integer,                             intent(in)    :: n, m
            type(surface_parcel_container_type), intent(inout) :: sp

            sp%position(n) = sp%position(m)

            sp%vorticity(n) = sp%vorticity(m)

            sp%length(n)  = sp%length(m)
            sp%buoyancy(n) = sp%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            sp%humidity(n) = sp%humidity(m)
#endif
        end subroutine surface_parcel_replace

        subroutine alloc(num, sp)
            integer,                             intent(in)    :: num
            type(surface_parcel_container_type), intent(inout) :: sp

            allocate(sp%position(num))
            allocate(sp%vorticity(num))
            allocate(sp%length(num))
            allocate(sp%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(sp%humidity(num))
#endif

            allocate(sp%delta_pos(num))
            allocate(sp%delta_vor(num))
            allocate(sp%strain(num))
        end subroutine alloc

        subroutine dealloc(sp)
            type(surface_parcel_container_type), intent(inout) :: sp

            deallocate(sp%position)
            deallocate(sp%vorticity)
            deallocate(sp%length)
            deallocate(sp%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(sp%humidity)
#endif
            deallocate(sp%delta_pos)
            deallocate(sp%delta_vor)
            deallocate(sp%strain)
        end subroutine dealloc

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine surface_parcel_alloc(num)
            integer, intent(in) :: num

            call alloc(num, top_parcels)
            call alloc(num, bot_parcels)

        end subroutine surface_parcel_alloc

        ! Deallocate parcel memory
        subroutine surface_parcel_dealloc

            call dealloc(top_parcels)
            call dealloc(bot_parcels)

        end subroutine surface_parcel_dealloc

end module surface_parcel_container
