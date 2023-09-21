! =============================================================================
! This module stores the surface parcel data and provides subroutines to write,
! modify allocate and deallocate it.
! =============================================================================
module surface_parcel_container
    use options, only : verbose
    use parameters, only : extent, extenti, center, lower, upper
    implicit none

    integer :: n_up_surf_parcels, n_lo_surf_parcels

    type surface_parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &   ! (x, y)
            vorticity,  &   ! (xi, eta, zeta)
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12, B(:, 3) = B22

        double precision, allocatable, dimension(:) :: &
            area,       &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy

            ! low-storage RK arrays:
        double precision, allocatable, dimension(:, :) :: &
            delta_pos,  &       ! velocity
            delta_vor,  &       ! vorticity tendency
            strain,     &
            delta_b             ! B-matrix tendency
    end type surface_parcel_container_type

    type(surface_parcel_container_type) up_surf_parcels, lo_surf_parcels

    private :: dealloc

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
            delx = delx - extent(1) * dble(nint(delx * extenti(1)))
        end function get_delx

#ifndef NDEBUG
        function get_dely(y1, y2) result (dely)
#else
        elemental function get_dely(y1, y2) result (dely)
#endif
            double precision, intent(in) :: y1, y2
            double precision             :: dely

            dely = y1 - y2
#ifndef NDEBUG
            if ((y1 < lower(2)) .or. (y2 < lower(2)) .or. (y1 > upper(2)) .or. (y2 > upper(2))) then
                write(*,*) 'point outside domain was fed into get_dely'
                write(*,*) 'y1, y2, lower(2), upper(2)'
                write(*,*) y1, y2, lower(2), upper(2)
            endif
#endif
            ! works across periodic edge
            dely = dely - extent(2) * dble(nint(dely * extenti(2)))
        end function get_dely


        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine surface_parcel_replace(s_parcels, n, m)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer, intent(in)                                :: n, m

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print '(a19, i0, a6, i0)', '    replace parcel ', n, ' with ', m
            endif
#endif

            s_parcels%position(:, n) = s_parcels%position(:, m)
            s_parcels%vorticity(:, n) = s_parcels%vorticity(:, m)

            s_parcels%area(n)  = s_parcels%area(m)
            s_parcels%B(:, n)  = s_parcels%B(:, m)

        end subroutine surface_parcel_replace

        ! Deallocate parcel memory
        subroutine surface_parcel_dealloc
            call dealloc(up_surf_parcels)
            call dealloc(lo_surf_parcels)
        end subroutine surface_parcel_dealloc


        subroutine surface_parcel_alloc(s_parcels, num)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(in)    :: num
            allocate(s_parcels%position(2, num))
            allocate(s_parcels%vorticity(3, num))
            allocate(s_parcels%B(3, num))
            allocate(s_parcels%area(num))
            allocate(s_parcels%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(s_parcels%humidity(num))
#endif
            ! low-storage RK arrays:
            allocate(s_parcels%delta_pos(2, num))
            allocate(s_parcels%delta_vor(3, num))
            allocate(s_parcels%strain(4, num))
            allocate(s_parcels%delta_b(3, num))
        end subroutine surface_parcel_alloc

        subroutine dealloc(s_parcels)
            type(surface_parcel_container_type), intent(inout) :: s_parcels

            deallocate(s_parcels%position)
            deallocate(s_parcels%vorticity)
            deallocate(s_parcels%B)
            deallocate(s_parcels%area)
            deallocate(s_parcels%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(s_parcels%humidity)
#endif
            ! low-storage RK arrays:
            deallocate(s_parcels%delta_pos)
            deallocate(s_parcels%delta_vor)
            deallocate(s_parcels%strain)
            deallocate(s_parcels%delta_b)
        end subroutine dealloc

end module surface_parcel_container
