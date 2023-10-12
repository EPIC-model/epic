! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_damping
    use constants, only : zero, one, two, f14
    use mpi_timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin, l_bndry_zeta_zero
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipsoid
    use parcel_interpl
    use fields
    use field_mpi, only : field_mpi_alloc                   &
                        , field_mpi_dealloc                 &
                        , field_buffer_to_halo              &
                        , field_halo_to_buffer              &
                        , field_buffer_to_interior          &
                        , field_interior_to_buffer          &
                        , interior_to_halo_communication    &
                        , halo_to_interior_communication    &
                        , field_halo_swap_scalar            
    use physics, only : glat, lambda_c, q_0
    use omp_lib
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)
    double precision :: weight(0:1,0:1,0:1)
    double precision :: time_fact(0:1,0:1,0:1)

    integer, parameter :: IDX_VOL_SWAP    = 1   &
                        , IDX_VORP_X_SWAP  = 2   &
                        , IDX_VORP_Y_SWAP  = 3   &
                        , IDX_VORP_Z_SWAP  = 4   &
                        , IDX_THETAP_SWAP  = 5   
#ifndef ENABLE_DRY_MODE
    integer, parameter :: IDX_QVP_SWAP   = 6
    integer, parameter :: IDX_QLP_SWAP   = 7

    integer, parameter :: n_field_swap = 7
#else
    integer, parameter :: n_field_swap = 5
#endif

    public :: parcel_damp

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt
            call par2grid_diag !get all gridded fields, including qv, ql and theta
            call perturbation_damping(dt)
        end subroutine parcel_damp

        ! 
        ! @pre 
        subroutine perturbation_damping(dt)
            double precision, intent(in)  :: dt
            integer                       :: n, p, l, ii, jj, kk
            double precision, parameter   :: pre_fact=1.0
            double precision              :: points(3, 4)
            double precision              :: pvol

            !call start_timer(grid2par_timer)

!            !$omp parallel default(shared)
!            !$omp do private(n, l, is, js, ks, weight, weights, reduce_fact, parcel_strain_mag) &
            do n = 1, n_parcels
                pvol = parcels%volume(n)
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n), n, .false.)

                ! we have 4 points per ellipsoid
                do p = 1, 4
                    call trilinear(points(:, p), is, js, ks, weights)
                    weight=0.25*weights
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
                    parcels%theta(n)=parcels%theta(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*thetag(ks:ks+1, js:js+1, is:is+1))
#ifndef ENABLE_DRY_MODE
                    parcels%qv(n)=parcels%qv(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*qvg(ks:ks+1, js:js+1, is:is+1))
                    parcels%ql(n)=parcels%ql(n)*(1.0-sum(weight*time_fact))&
                    +sum(weight*time_fact*qlg(ks:ks+1, js:js+1, is:is+1))
                enddo
#endif
             enddo

        end subroutine perturbation_damping


end module parcel_damping
