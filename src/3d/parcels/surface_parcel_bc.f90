! =============================================================================
!                       Surface parcel boundary conditions
!                       periodic in x and y (horizontal)
! =============================================================================
module surface_parcel_bc
    use parameters, only : extent, hli, center
    use surface_parcel_container, only : lo_surf_parcels, n_lo_surf_parcels &
                                       , up_surf_parcels, n_up_surf_parcels
    use omp_lib
    implicit none

    private :: apply_parcel_bc

    contains

        ! Apply periodic bc on n-th parcel (horizontally)
        ! @param[inout] position vector of parcel
        subroutine apply_surface_periodic_bc(position)
            double precision, intent(inout) :: position(2)
            position(1) = position(1) - extent(1) * dble(int((position(1) - center(1)) * hli(1)))
            position(2) = position(2) - extent(2) * dble(int((position(2) - center(2)) * hli(2)))
        end subroutine apply_surface_periodic_bc

        subroutine apply_surface_parcel_bc
            call apply_parcel_bc(lo_surf_parcels%position, n_lo_surf_parcels)
            call apply_parcel_bc(up_surf_parcels%position, n_up_surf_parcels)
        end subroutine apply_surface_parcel_bc

        ! Apply all boundary conditions to all parcels
        ! @param[inout] position vector of parcels
        ! @param[inout] B matrix of parcels
        subroutine apply_parcel_bc(position, n_par)
            double precision, intent(inout) :: position(:, :)
            integer,          intent(in)    :: n_par
            integer                         :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                ! horizontal direction
                call apply_surface_periodic_bc(position(:, n))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine apply_parcel_bc

end module surface_parcel_bc
