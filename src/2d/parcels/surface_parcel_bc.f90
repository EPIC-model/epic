! =============================================================================
!                       Survface parcel boundary conditions
!                       periodic in x (horizontal)
! =============================================================================
module surface_parcel_bc
    use parameters, only : extent, hli, center
    use surface_parcel_container, only : n_surf_parcels
    use omp_lib
    implicit none

    contains

        ! Apply periodic bc on n-th parcel (horizontally)
        ! @param[inout] position vector of parcel
        subroutine apply_surface_periodic_bc(position)
            double precision, intent(inout) :: position
            position = position - extent(1) * dble(int((position - center(1)) * hli(1)))
        end subroutine apply_surface_periodic_bc

        ! Apply all boundary conditions to all parcels
        ! @param[inout] position vector of parcels
        subroutine apply_surface_parcel_bc(position)
            double precision, intent(inout) :: position(:)
            integer                         :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_surf_parcels
                ! horizontal direction
                call apply_surface_periodic_bc(position(n))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine apply_surface_parcel_bc

end module surface_parcel_bc
