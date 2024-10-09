! =============================================================================
!                       Parcel boundary conditions
!                       periodic in x (zonal) and in y (meridional)
!                       reflective in z (vertical)
! =============================================================================
module parcel_bc
    use constants, only : zero, two
    use parameters, only : lower, upper, extent, extenti
    use parcel_container, only : n_parcels, parcels
    use omp_lib
    implicit none

    contains

        ! Apply periodic bc on n-th parcel (zonally and meridionally)
        ! @param[inout] position vector of parcel
        pure subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(3)

            if(position(1) >= upper(1)) then
                position(1) = position(1) - int(1+(position(1)-upper(1))*extenti(1))*extent(1)
            else if (position(1) < lower(1)) then
                position(1) = position(1) + int(1+(lower(1)-position(1))*extenti(1))*extent(1)
            endif
            
            if(position(2) >= upper(2)) then
                position(2) = position(2) - int(1+(position(2)-upper(2))*extenti(2))*extent(2)
            else if (position(2) < lower(2)) then
                position(2) = position(2) + int(1+(lower(2)-position(2))*extenti(2))*extent(2)
            endif

        end subroutine apply_periodic_bc

        ! Apply mirroring bc on n-th parcel (vertically)
        ! @param[inout] position vector of parcel
        ! @param[inout] B matrix of parcel
        pure subroutine apply_reflective_bc(position, B)
            double precision, intent(inout) :: position(3), B(5)

            if (position(3) > upper(3)) then
                position(3) = two * upper(3) - position(3)
                ! flip sign of B13 and B23
                B(3) = -B(3)
                B(5) = -B(5)
            else if (position(3) < lower(3)) then
                position(3) = two * lower(3) - position(3)
                ! flip sign of B13 and B23
                B(3) = -B(3)
                B(5) = -B(5)
            endif
        end subroutine apply_reflective_bc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Apply reflective boundary conditions to all parcels
        subroutine apply_parcel_reflective_bc
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! vertical direction
                call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_parcel_reflective_bc

end module parcel_bc
