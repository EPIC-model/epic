! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_damping
    use constants, only :  f14
    use mpi_timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin 
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid
    use parcel_interpl
    use fields
    use omp_lib
    implicit none

    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)
    double precision :: weight(0:1,0:1,0:1)
    double precision :: time_fact(0:1,0:1,0:1)

    public :: parcel_damp

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt
            call par2grid 
            call perturbation_damping(dt)
        end subroutine parcel_damp

        ! 
        ! @pre 
        subroutine perturbation_damping(dt)
            double precision, intent(in)  :: dt
            integer                       :: n, p, l, ii, jj, kk
            double precision, parameter   :: pre_fact=1.0
            double precision              :: points(3,4)
            double precision              :: pvol

            !call start_timer(grid2par_timer)

!           !$omp parallel default(shared)
!           !$omp do private(n, l, is, js, ks, weight, weights, time_fact, pvol, points)&
!           !$omp do private&(p, ii, jj ,kk) 
            do n = 1, n_parcels
                pvol = parcels%volume(n)
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n), n, .false.)

                ! we have 4 points per ellipsoid
                do p = 1, 4
                    call trilinear(points(:, p), is, js, ks, weights)
                    weight=f14*weights
                    do ii=0,1
                      do jj=0,1
                        do kk=0,1
                           time_fact(kk,jj,ii)=1.0-exp(-pre_fact*strain_mag(ks+kk, js+jj, is+ii)*dt)
                        enddo
                      enddo
                    enddo
                    do l=1,3
                         parcels%vorticity(l,n)=parcels%vorticity(l,n)*(1.0-sum(weight*time_fact))&
                         +sum(weight*time_fact*vortg(ks:ks+1, js:js+1, is:is+1, l))
                    enddo
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n)=parcels%humidity(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*humg(ks:ks+1, js:js+1, is:is+1))
                    parcels%buoyancy(n)=parcels%buoyancy(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*dbuoyg(ks:ks+1, js:js+1, is:is+1))
#else
                    parcels%buoyancy(n)=parcels%buoyancy(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*tbuoyg(ks:ks+1, js:js+1, is:is+1))
#endif
                enddo
            enddo

        end subroutine perturbation_damping


end module parcel_damping
