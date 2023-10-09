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
                        , field_halo_swap_scalar            &
                        , field_halo_to_buffer_integer      &
                        , field_buffer_to_interior_integer  &
                        , field_interior_to_buffer_integer  &
                        , field_buffer_to_halo_integer
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

    integer, parameter :: IDX_VOL_SWAP    = 1   &
                        , IDX_VORP_X_SWAP  = 2   &
                        , IDX_VORP_Y_SWAP  = 3   &
                        , IDX_VORP_Z_SWAP  = 4   &
                        , IDX_BUOYP_SWAP  = 5   
#ifndef ENABLE_DRY_MODE
    integer, parameter :: IDX_HUMP_SWAP   = 6

    integer, parameter :: n_field_swap = 6
#else
    integer, parameter :: n_field_swap = 5
#endif

    public :: parcel_damp

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt
            call par2grid(.false.)
            call perturbation_damping(dt)
        end subroutine parcel_damp

        ! 
        ! @pre 
        subroutine perturbation_damping(dt)
            double precision, intent(in)  :: dt
            integer                       :: n, l
            double precision              :: parcel_strain_mag
            double precision              :: parcel_vortp(3)
            double precision              :: parcel_hump
            double precision              :: parcel_buoyp
            double precision              :: reduce_fact
            double precision, parameter   :: pre_fact=0.5

            !call start_timer(grid2par_timer)

            vortpg=zero
#ifndef ENABLE_DRY_MODE
            humpg=zero
#endif
            buoypg=zero
            volg=zero

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, ks, weight, weights, reduce_fact) &
#ifndef ENABLE_DRY_MODE
            !$omp& private(parcel_strain_mag, parcel_vortp, parcel_hump, parcel_buoyp)& 
            !$omp& reduction(+: vortpg, humpg, buoypg, volg)
#else
            !$omp& private(parcel_strain_mag, parcel_vortp, parcel_buoyp)& 
            !$omp& reduction(+: vortpg, buoypg, volg)
#endif

            do n = 1, n_parcels

                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(:, n), is, js, ks, weights)

                ! loop over grid points which are part of the interpolation
                parcel_strain_mag=sum(weights * strain_mag(ks:ks+1, js:js+1, is:is+1))
                reduce_fact=(1.0-exp(-pre_fact*parcel_strain_mag*dt))
                do l=1,3
                    parcel_vortp(l)=reduce_fact*parcels%vorticity(l, n)
                enddo
#ifndef ENABLE_DRY_MODE
                parcel_hump=reduce_fact*parcels%humidity(n)
                parcel_buoyp=reduce_fact*parcels%buoyancy(n)
#else
                parcel_buoyp=reduce_fact*parcels%buoyancy(n)
#endif
                weight = parcels%volume(n) * weights

                do l = 1, 3
                    vortpg(ks:ks+1, js:js+1, is:is+1, l) = vortpg(ks:ks+1, js:js+1, is:is+1, l) &
                                                            + weight *parcel_vortp(l)
                enddo
                buoypg(ks:ks+1, js:js+1, is:is+1) = buoypg(ks:ks+1, js:js+1, is:is+1) &
                                                      + weight * parcel_buoyp
#ifndef ENABLE_DRY_MODE
                humpg(ks:ks+1, js:js+1, is:is+1) = humpg(ks:ks+1, js:js+1, is:is+1) &
                                              + weight * parcel_hump
#endif
                volg(ks:ks+1, js:js+1, is:is+1) = volg(ks:ks+1, js:js+1, is:is+1) &
                                              + weight *reduce_fact
                 
            enddo
            !$omp end do
            !$omp end parallel

            ! DO THE HALO STUFF
            call perturbation_damping_halo_swap 

            !$omp parallel workshare
            ! calculate the average correction needed at grid point
            ! exclude halo cells to avoid division by zero
            do l = 1, 3
                vortpg(0:nz, :, :, l) = vortpg(0:nz, :, :, l) / volg(0:nz, :, :)
            enddo


#ifndef ENABLE_DRY_MODE
            humpg(0:nz, :, :) = humpg(0:nz, :, :) / volg(0:nz, :, :)
#endif
            buoypg(0:nz, :, :) = buoypg(0:nz, :, :) / volg(0:nz, :, :)
            !$omp end parallel workshare

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, ks, weights, parcel_strain_mag, reduce_fact)
            do n = 1, n_parcels
                call trilinear(parcels%position(:, n), is, js, ks, weights)
                parcel_strain_mag=sum(weights * strain_mag(ks:ks+1, js:js+1, is:is+1))
                reduce_fact=(1.0-exp(-pre_fact*parcel_strain_mag*dt))
                do l=1,3
                    parcels%vorticity(l,n)=parcels%vorticity(l,n)*(1.0-reduce_fact)+reduce_fact*&
                                           sum(weights * vortpg(ks:ks+1, js:js+1, is:is+1, l))
                enddo
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n)=parcels%humidity(n)*(1.0-reduce_fact)+reduce_fact*&
                                    sum(weights * humpg(ks:ks+1, js:js+1, is:is+1))
#endif
                parcels%buoyancy(n)=parcels%buoyancy(n)*(1.0-reduce_fact)+reduce_fact*&
                                    sum(weights * buoypg(ks:ks+1, js:js+1, is:is+1))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine perturbation_damping


        subroutine perturbation_damping_halo_swap
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_mpi_alloc(n_field_swap, ndim=3)

            !------------------------------------------------------------------
            ! Accumulate interior:

            call field_halo_to_buffer(volg,                 IDX_VOL_SWAP)
            call field_halo_to_buffer(vortpg(:, :, :, I_X), IDX_VORP_X_SWAP)
            call field_halo_to_buffer(vortpg(:, :, :, I_Y), IDX_VORP_Y_SWAP)
            call field_halo_to_buffer(vortpg(:, :, :, I_Z), IDX_VORP_Z_SWAP)
            call field_halo_to_buffer(buoypg,               IDX_BUOYP_SWAP)
#ifndef ENABLE_DRY_MODE
            call field_halo_to_buffer(humpg,                IDX_HUMP_SWAP)
#endif

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call field_buffer_to_interior(volg,                 IDX_VOL_SWAP, .true.)
            call field_buffer_to_interior(vortpg(:, :, :, I_X), IDX_VORP_X_SWAP, .true.)
            call field_buffer_to_interior(vortpg(:, :, :, I_Y), IDX_VORP_Y_SWAP, .true.)
            call field_buffer_to_interior(vortpg(:, :, :, I_Z), IDX_VORP_Z_SWAP, .true.)
            call field_buffer_to_interior(buoypg,               IDX_BUOYP_SWAP, .true.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_interior(humpg,                IDX_HUMP_SWAP, .true.)
#endif

            !------------------------------------------------------------------
            ! Fill halo:

            call field_interior_to_buffer(volg,                 IDX_VOL_SWAP)
            call field_interior_to_buffer(vortpg(:, :, :, I_X), IDX_VORP_X_SWAP)
            call field_interior_to_buffer(vortpg(:, :, :, I_Y), IDX_VORP_Y_SWAP)
            call field_interior_to_buffer(vortpg(:, :, :, I_Z), IDX_VORP_Z_SWAP)
            call field_interior_to_buffer(buoypg,               IDX_BUOYP_SWAP)
#ifndef ENABLE_DRY_MODE
            call field_interior_to_buffer(humpg,                IDX_HUMP_SWAP)
#endif

            call interior_to_halo_communication

            call field_buffer_to_halo(volg,                 IDX_VOL_SWAP, .false.)
            call field_buffer_to_halo(vortpg(:, :, :, I_X), IDX_VORP_X_SWAP, .false.)
            call field_buffer_to_halo(vortpg(:, :, :, I_Y), IDX_VORP_Y_SWAP, .false.)
            call field_buffer_to_halo(vortpg(:, :, :, I_Z), IDX_VORP_Z_SWAP, .false.)
            call field_buffer_to_halo(buoypg,               IDX_BUOYP_SWAP, .false.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_halo(humpg,                 IDX_HUMP_SWAP, .false.)
#endif

            call field_mpi_dealloc

        end subroutine perturbation_damping_halo_swap


end module parcel_damping
