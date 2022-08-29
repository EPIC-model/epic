! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, hli, center, lower, upper
    use parcel_ellipsoid, only : parcel_ellipsoid_allocate, parcel_ellipsoid_deallocate
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

    integer, parameter :: I_B11 = 1 & ! index for B11 matrix component
                        , I_B12 = 2 & ! index for B12 matrix component
                        , I_B13 = 3 & ! index for B13 matrix component
                        , I_B22 = 4 & ! index for B22 matrix component
                        , I_B23 = 5   ! index for B23 matrix component

    contains

        ! Obtain the difference between two zonal coordinates
        ! across periodic edges
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
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

        ! Obtain the difference between two meridional coordinates
        ! across periodic edges
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        ! WARNING input needs to be between lower and upper (see debug statement)
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
            dely = dely - extent(2) * dble(int(dely * hli(2)))
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
            call parcel_ellipsoid_allocate
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
            call parcel_ellipsoid_deallocate
        end subroutine parcel_dealloc

end module parcel_container
