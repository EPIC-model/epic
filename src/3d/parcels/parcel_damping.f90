! =============================================================================
! This module contains the subroutines to do damping of parcel properties to
! gridded fields in a conservative manner.
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
    use options, only : damping
    implicit none

    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks
    integer :: damping_timer
 
    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)
    double precision :: weight(0:1,0:1,0:1)
    double precision :: time_fact(0:1,0:1,0:1)

    public :: parcel_damp, damping_timer

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt

            if(damping%l_vorticity .or. damping%l_scalars) then 
                call par2grid(.false.)
                call perturbation_damping(dt,damping%l_vorticity,damping%l_scalars,&
                     damping%vorticity_prefactor,damping%scalars_prefactor,.true.)
            end if

        end subroutine parcel_damp

        ! 
        ! @pre 
        subroutine perturbation_damping(dt,l_vorticity,l_scalars,&
                                        vorticity_prefactor,scalars_prefactor,l_reuse)
            double precision, intent(in)  :: dt
            logical, intent(in)   :: l_vorticity
            logical, intent(in)   :: l_scalars
            logical   :: l_reuse
            double precision, intent(in)   :: vorticity_prefactor
            double precision, intent(in)   :: scalars_prefactor
            integer                       :: n, p, l, ii, jj, kk
            double precision              :: points(3, n_points_p2g)
            double precision              :: pvol

            call start_timer(damping_timer)

            !$omp parallel default(shared)
            !$omp do private(n, p, l, ii, jj, kk, points, pvol, weight) &
            !$omp& private( is, js, ks, weights, time_fact) 
            do n = 1, n_parcels
                pvol = parcels%volume(n)
#ifndef ENABLE_P2G_1POINT
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n), n, l_reuse)
#else
                points(:, 1) = parcels%position(:, n)
#endif

                ! we have 4 points per ellipsoid
                do p = 1, n_points_p2g
                    call trilinear(points(:, p), is, js, ks, weights)
                    weight = point_weight_p2g * weights
                    if(l_vorticity) then
                      do ii=0,1
                        do jj=0,1
                          do kk=0,1
                            time_fact(kk,jj,ii)=1.0-exp(-vorticity_prefactor*&
                                                strain_mag(ks+kk, js+jj, is+ii)*dt)
                          enddo
                        enddo
                      enddo
                      do l=1,3
                         parcels%vorticity(l,n)=parcels%vorticity(l,n)*(1.0-sum(weight*time_fact))&
                         +sum(weight*time_fact*vortg(ks:ks+1, js:js+1, is:is+1, l))
                      enddo
                    endif
                    if(l_scalars) then
                      do ii=0,1
                        do jj=0,1
                          do kk=0,1
                            time_fact(kk,jj,ii)=1.0-exp(-scalars_prefactor*&
                                                strain_mag(ks+kk, js+jj,is+ii)*dt)
                          enddo
                        enddo
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
                    endif
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(damping_timer)


        end subroutine perturbation_damping


end module parcel_damping
