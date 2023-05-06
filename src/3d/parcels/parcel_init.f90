! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel
    use constants, only : zero, two, one, f12, f13, f23, f14
    use parcel_container, only : parcels, n_parcels, n_total_parcels, parcel_alloc
    use parcel_ellipsoid, only : get_abc, get_eigenvalues
    use parcel_split_mod, only : parcel_split
    use parcel_interpl, only : trilinear, ngp
    use parameters, only : dx, vcell, ncell,            &
                           extent, lower, nx, ny, nz,   &
                           max_num_parcels
    use mpi_timer, only : start_timer, stop_timer
    use field_mpi, only : field_halo_fill
    use omp_lib
    use mpi_communicator
    use mpi_layout, only : box
    use mpi_utils, only : mpi_print, mpi_exit_on_error
    use mpi_collectives, only : mpi_blocking_reduce
    use fields
    implicit none

    integer :: init_timer

    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(ngp)

    private :: weights, is, js, ks

    private :: init_refine

    contains

        ! Allocate parcel container and sets values for parcel attributes
        ! to their default values.
        subroutine parcel_default
            double precision             :: lam, l23
            integer                      :: n

            call start_timer(init_timer)

            call parcel_alloc(max_num_parcels)

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel%n_per_cell * box%ncell

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                call mpi_exit_on_error
            endif

            n_total_parcels = n_parcels
            if (comm%size > 1) then
                call mpi_blocking_reduce(n_total_parcels, MPI_SUM)
            endif

            call init_regular_positions

            ! initialize the volume of each parcel
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%volume(n) = vcell / dble(parcel%n_per_cell)
            enddo
            !$omp end do
            !$omp end parallel

            ! aspect ratio: lam = a / c
            lam = maxval((/dx(1) / dx(2), dx(2) / dx(1),   &
                           dx(1) / dx(3), dx(3) / dx(1),   &
                           dx(2) / dx(3), dx(3) / dx(2)/))

            !$omp parallel default(shared)
            !$omp do private(n, l23)
            do n = 1, n_parcels
                ! set all to zero
                parcels%B(:, n) = zero

                l23 = (lam * get_abc(parcels%volume(n))) ** f23

                ! B11
                parcels%B(1, n) = l23

                ! B22
                parcels%B(4, n) = l23
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(lam)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%vorticity(:, n) = zero
                parcels%buoyancy(n) = zero
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n) = zero
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(init_timer)

        end subroutine parcel_default


        ! Position parcels regularly in the domain.
        subroutine init_regular_positions
            integer          :: ix, i, iz, j, iy, k, l, n_per_dim
            double precision :: im, corner(3)

            ! number of parcels per dimension
            n_per_dim = int(dble(parcel%n_per_cell) ** f13)
            if (n_per_dim ** 3 .ne. parcel%n_per_cell) then
                if (comm%rank == comm%master) then
                    print *, "Number of parcels per cell (", &
                             parcel%n_per_cell, ") not a cubic."
                endif
                call mpi_exit_on_error
            endif

            im = one / dble(n_per_dim)

            l = 1
            do iz = 0, nz-1
                do iy = box%lo(2), box%hi(2)
                    do ix = box%lo(1), box%hi(1)
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do k = 1, n_per_dim
                            do j = 1, n_per_dim
                                do i = 1, n_per_dim
                                    parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                                    parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                                    parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                                    l = l + 1
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_parcels == l - 1) then
                call mpi_exit_on_error("Number of parcels disagree!")
            endif
        end subroutine init_regular_positions


        subroutine init_refine(lam)
            double precision, intent(inout) :: lam
            double precision                :: evals(3) ! = (a2, b2, c2)

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call parcel_split(parcels, parcel%lambda_max)
                evals = get_eigenvalues(parcels%B(1, :), parcels%volume(1))
                lam = dsqrt(evals(1) / evals(3))
            end do
        end subroutine init_refine


        subroutine init_parcels_from_grids
            integer:: n, l

            call start_timer(init_timer)

            ! make usre halo grid points are filled
            call field_halo_fill(vortg(:, :, :, I_X))
            call field_halo_fill(vortg(:, :, :, I_Y))
            call field_halo_fill(vortg(:, :, :, I_Z))
            call field_halo_fill(tbuoyg)
#ifndef ENABLE_DRY_MODE
            call field_halo_fill(humg)
#endif

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, ks, weights)
            do n = 1, n_parcels

                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(:, n), is, js, ks, weights)

                ! loop over grid points which are part of the interpolation
                parcels%vorticity(:, n) = parcels%vorticity(:, n) &
                                   + weights(1) * vortg(ks, js, is, :) &
                                   + weights(2) * vortg(ks  , js   , is+1, :) &
                                   + weights(3) * vortg(ks  , js+1 , is, :) &
                                   + weights(4) * vortg(ks  , js+1 , is+1, :) &
                                   + weights(5) * vortg(ks+1, js   , is, :) &
                                   + weights(6) * vortg(ks+1, js   , is+1, :) &
                                   + weights(7) * vortg(ks+1, js+1 , is, :) &
                                   + weights(8) * vortg(ks+1, js+1 , is+1, :)
                parcels%buoyancy(n) = parcels%buoyancy(n) &
                                   + weights(1) * tbuoyg(ks, js, is) &
                                   + weights(2) * tbuoyg(ks  , js   , is+1) &
                                   + weights(3) * tbuoyg(ks  , js+1 , is) &
                                   + weights(4) * tbuoyg(ks  , js+1 , is+1) &
                                   + weights(5) * tbuoyg(ks+1, js   , is) &
                                   + weights(6) * tbuoyg(ks+1, js   , is+1) &
                                   + weights(7) * tbuoyg(ks+1, js+1 , is) &
                                   + weights(8) * tbuoyg(ks+1, js+1 , is+1)
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n) = parcels%humidity(n) &
                                   + weights(1) * humg(ks, js, is) &
                                   + weights(2) * humg(ks  , js   , is+1) &
                                   + weights(3) * humg(ks  , js+1 , is) &
                                   + weights(4) * humg(ks  , js+1 , is+1) &
                                   + weights(5) * humg(ks+1, js   , is) &
                                   + weights(6) * humg(ks+1, js   , is+1) &
                                   + weights(7) * humg(ks+1, js+1 , is) &
                                   + weights(8) * humg(ks+1, js+1 , is+1)
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(init_timer)

        end subroutine init_parcels_from_grids

end module parcel_init
