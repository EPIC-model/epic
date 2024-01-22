! =============================================================================
!                       Survface parcel boundary conditions
!                       periodic in x (horizontal)
! =============================================================================
module surface_parcel_bc
    use parameters, only : extent, hli, center
    use omp_lib
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_top_parcels, top_parcels     &
                                       , n_bot_parcels, bot_parcels
    implicit none

    private :: apply_surface_parcel_bc_

    contains

        ! Apply periodic bc on n-th parcel (horizontally)
        ! @param[inout] position vector of parcel
        subroutine apply_surface_periodic_bc(position)
            double precision, intent(inout) :: position
            position = position - extent(1) * dble(int((position - center(1)) * hli(1)))
        end subroutine apply_surface_periodic_bc

        subroutine apply_surface_parcel_bc

            call apply_surface_parcel_bc_(n_bot_parcels, bot_parcels)
            call apply_surface_parcel_bc_(n_top_parcels, top_parcels)

        end subroutine apply_surface_parcel_bc

        ! Apply all boundary conditions to all parcels
        ! @param[inout] position vector of parcels
        subroutine apply_surface_parcel_bc_(n_par, spar)
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer                                            :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                ! horizontal direction
                call apply_surface_periodic_bc(spar%position(n))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine apply_surface_parcel_bc_

end module surface_parcel_bc
