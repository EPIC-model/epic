!===================================================================
!       Module with divergent flow and gradient correction
!===================================================================
module parcel_correction
    use inversion_utils, only : init_inversion
    use inversion_mod, only : diverge
    use parcel_interpl, only : trilinear, ngp, vol2grid
    use parcel_bc
    use omp_lib
    use constants
    use parameters, only : vcelli, nx, ny, nz, dx
    use parcel_container
    use mpi_timer, only : start_timer, stop_timer
    use fields, only : volg
    use mpi_layout, only : box
    use mpi_communicator
    use mpi_utils, only : mpi_check_for_error
    implicit none

    private

    integer :: lapl_corr_timer, &
               grad_corr_timer, &
               vort_corr_timer

    ! initial volume-weighted vorticity mean
    double precision :: vor_bar(3)

    public :: init_parcel_correction, &
              apply_laplace,          &
              apply_gradient,         &
              apply_vortcor,          &
              lapl_corr_timer,        &
              grad_corr_timer,        &
              vort_corr_timer

    contains

        ! Initialise parcel correction module
        subroutine init_parcel_correction
            integer          :: n
            double precision :: vsum
            double precision :: buf(4)

            call start_timer(vort_corr_timer)

            call init_inversion

            vsum = zero
            vor_bar = zero

            !$omp parallel default(shared)
            !$omp do private(n) reduction(+: vor_bar, vsum)
            do n = 1, n_parcels
                vsum = vsum + parcels%volume(n)
                vor_bar = vor_bar + parcels%vorticity(:, n) * parcels%volume(n)
            enddo
            !$omp end do
            !$omp end parallel

            buf(1) = vsum
            buf(2:4) = vor_bar
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:4),                &
                               4,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               comm%world,              &
                               comm%err)

            call mpi_check_for_error("in MPI_Allreduce of parcel_correction::init_parcel_correction.")

            vsum = buf(1)
            vor_bar = buf(2:4)

            vor_bar = vor_bar / vsum

            call stop_timer(vort_corr_timer)

        end subroutine init_parcel_correction

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This subroutine makes sure that the net vorticity
        ! (volume-integrated vorticity) remains constant.
        subroutine apply_vortcor
            integer          :: n
            double precision :: dvor(3), vsum
            double precision :: buf(4)

            call start_timer(vort_corr_timer)

            vsum = zero
            dvor = - vor_bar

            !$omp parallel default(shared)
            !$omp do private(n) reduction(+: dvor, vsum)
            do n = 1, n_parcels
                vsum = vsum + parcels%volume(n)
                dvor = dvor + parcels%vorticity(:, n) * parcels%volume(n)
            enddo
            !$omp end do
            !$omp end parallel

            buf(1) = vsum
            buf(2:4) = dvor
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:4),                &
                               4,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               comm%world,              &
                               comm%err)

            call mpi_check_for_error("in MPI_Allreduce of parcel_correction::apply_vortcor.")
            vsum = buf(1)
            dvor = buf(2:4)

            dvor = dvor / vsum

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%vorticity(:, n) = parcels%vorticity(:, n) - dvor
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(vort_corr_timer)
        end subroutine apply_vortcor

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_laplace(l_reuse)
            logical, optional, intent(in) :: l_reuse
            double precision              :: phi(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)), &
                                              ud(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)), &
                                              vd(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)), &
                                              wd(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision              :: weights(0:1,0:1,0:1)
            integer                       :: n, l, is, js, ks

            call start_timer(lapl_corr_timer)

            call vol2grid(l_reuse)

            ! form divergence field
            phi = volg * vcelli - one

            call diverge(phi, ud, vd, wd)

            ! use extrapolation to fill z grid lines outside domain:
            ! (anti-symmetry in w is equal to extrapolation)
            ud(-1, :, :)   =  two * ud(0, :, :) - ud(1, :, :)
            vd(-1, :, :)   =  two * vd(0, :, :) - vd(1, :, :)
            wd(-1, :, :)   = -wd(1, :, :)
            ud(nz+1, :, :) =  two * ud(nz, :, :) - ud(nz-1, :, :)
            vd(nz+1, :, :) =  two * vd(nz, :, :) - vd(nz-1, :, :)
            wd(nz+1, :, :) = -wd(nz-1, :, :)

            !------------------------------------------------------------------
            ! Increment parcel positions usind (ud, vd, wd) field:
            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights)
            do n = 1, n_parcels
                call trilinear(parcels%position(:, n), is, js, ks, weights)

                parcels%position(1, n) = parcels%position(1, n)               &
                                           + sum(weights(:,:,:) * ud(ks:ks+1, js:js+1, is:is+1))

                parcels%position(2, n) = parcels%position(2, n)               &
                                           + sum(weights(:,:,:) * vd(ks:ks+1, js:js+1, is:is+1))

                parcels%position(3, n) = parcels%position(3, n)               &
                                           + sum(weights(:,:,:) * wd(ks:ks+1, js:js+1, is:is+1))
            enddo
            !$omp end do
            !$omp end parallel

            call apply_swap_periodicity

            call stop_timer(lapl_corr_timer)

        end subroutine apply_laplace

        subroutine apply_gradient(prefactor, max_compression, l_reuse)
            double precision,  intent(in) :: prefactor
            double precision,  intent(in) :: max_compression
            logical, optional, intent(in) :: l_reuse
            double precision              :: phi(0:nz, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision              :: weights(0:1,0:1,0:1)
            double precision              :: xs, ys, zs, xf, yf, zf, xfc, yfc, zfc, lim_x, lim_y, lim_z
            integer                       :: n, is, js, ks

            call start_timer(grad_corr_timer)

            call vol2grid(l_reuse)

            ! form divergence field * dt and store in phi temporarily:
            phi = volg(0:nz, :, :) * vcelli - one

            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights, xf, yf, zf, xfc, yfc, zfc, xs, ys, zs, lim_x, lim_y, lim_z)
            do n = 1, n_parcels

                call trilinear(parcels%position(:, n), is, js, ks, weights)

                xf = weights(0,0,1) + weights(0,1,1) + weights(1,0,1) + weights(1,1,1) ! fractional position along x
                yf = weights(0,1,0) + weights(0,1,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along y
                zf = weights(1,0,0) + weights(1,0,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along z
                xfc=one-xf
                yfc=one-yf
                zfc=one-zf

                ! l   ks(l) js(l) is(l)
                ! 1   k     j     i
                ! 2   k     j     i+1
                ! 3   k     j+1   i
                ! 4   k     j+1   i+1
                ! 5   k+1   j     i
                ! 6   k+1   j     i+1
                ! 7   k+1   j+1   i
                ! 8   k+1   j+1   i+1
                lim_x = dx(1) * xf * xfc
                xs = - prefactor * lim_x * (&
                            zfc * (&
                                    yfc * (phi(ks, js,   is+1) - phi(ks, js,   is))  &
                                  +       yf * (phi(ks, js+1, is+1) - phi(ks, js+1, is))) &
                          +      zf * (&
                                    yfc * (phi(ks+1, js,   is+1) - phi(ks+1, js,   is))  &
                                  +       yf * (phi(ks+1, js+1, is+1) - phi(ks+1, js+1, is))))

                lim_x = lim_x * max_compression
                xs = max(-lim_x, min(xs, lim_x))

                lim_y = dx(2) * yf * yfc
                ys = - prefactor * lim_y * (&
                            zfc * (&
                                    xfc * (phi(ks, js+1, is  ) - phi(ks, js, is))  &
                                  +       xf * (phi(ks, js+1, is+1) - phi(ks, js, is+1))) &
                          +      zf * (&
                                    xfc * (phi(ks+1, js+1, is) -   phi(ks+1, js, is))  &
                                  +       xf * (phi(ks+1, js+1, is+1) - phi(ks+1, js, is+1))))

                lim_y = lim_y * max_compression
                ys = max(-lim_y, min(ys, lim_y))

                lim_z = dx(3) * zf * zfc
                zs = - prefactor * lim_z * (&
                            xfc * (&
                                    yfc * (phi(ks+1, js, is) - phi(ks, js, is))  &
                                  +       yf * (phi(ks+1, js+1, is)) - phi(ks, js+1, is)) &
                          +      xf * (&
                                    yfc * (phi(ks+1, js, is+1) - phi(ks, js, is+1))  &
                                  +       yf * (phi(ks+1, js+1, is+1) - phi(ks, js+1, is+1))))

                lim_z = lim_z * max_compression
                zs = max(-lim_z, min(zs, lim_z))

                parcels%position(1, n) = parcels%position(1, n) + xs
                parcels%position(2, n) = parcels%position(2, n) + ys
                parcels%position(3, n) = parcels%position(3, n) + zs
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(grad_corr_timer)

        end subroutine apply_gradient

end module parcel_correction
