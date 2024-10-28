!===================================================================
!       Module with divergent flow and gradient correction
!===================================================================
module parcel_correction
    use inversion_utils, only : init_inversion
    use inversion_mod, only : diverge
    use parcel_interpl, only : vol2grid
    use interpl, only : trilinear
    use parcel_bc
    use omp_lib
    use constants
    use parameters, only : vcelli, nx, ny, nz, dx
    use parcels_mod, only : parcels
    use mpi_timer, only : start_timer, stop_timer
    use fields, only : volg
    use mpi_layout, only : box
    use mpi_environment
    use mpi_utils, only : mpi_check_for_error
    use parcel_mpi, only : parcel_communicate
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
            do n = 1, parcels%local_num
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
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of parcel_correction::init_parcel_correction.")

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
            dvor = zero

            !$omp parallel default(shared)
            !$omp do private(n) reduction(+: dvor, vsum)
            do n = 1, parcels%local_num
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
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of parcel_correction::apply_vortcor.")
            vsum = buf(1)
            dvor = buf(2:4)

            ! vorticity correction:
            ! vor_mean = dvor / vsum
            ! if vor_mean < vor_bar --> correction < 0 --> we must add to parcel vorticity
            ! if vor_mean > vor_bar --> correction > 0 --> we must subtract from parcel vorticity
            dvor = dvor / vsum - vor_bar

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, parcels%local_num
                parcels%vorticity(:, n) = parcels%vorticity(:, n) - dvor
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(vort_corr_timer)
        end subroutine apply_vortcor

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_laplace(l_reuse)
            logical, optional, intent(in) :: l_reuse
            double precision              :: phi(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision              :: grad(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 3)
            double precision              :: weights(0:1,0:1,0:1)
            integer                       :: n, l, is, js, ks

            call start_timer(lapl_corr_timer)

            call vol2grid(l_reuse)

            ! form divergence field
            phi = volg * vcelli - one

            call diverge(phi, grad)

            ! use extrapolation to fill z grid lines outside domain:
            grad(-1,   :, :, :) =  two * grad(0,  :, :, :) - grad(1,    :, :, :)
            grad(nz+1, :, :, :) =  two * grad(nz, :, :, :) - grad(nz-1, :, :, :)

            !------------------------------------------------------------------
            ! Increment parcel positions usind (ud, vd, wd) field:
            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, ks, weights)
            do n = 1, parcels%local_num
                call trilinear(parcels%position(:, n), is, js, ks, weights)

                do l = 1, 3
                    parcels%position(l, n) = parcels%position(l, n)                     &
                                           + sum(weights * grad(ks:ks+1, js:js+1, is:is+1, l))
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call parcel_communicate(parcels)

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
            phi = volg * vcelli - one

            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights, xf, yf, zf, xfc, yfc, zfc, xs, ys, zs, lim_x, lim_y, lim_z)
            do n = 1, parcels%local_num

                call trilinear(parcels%position(:, n), is, js, ks, weights)

                xf = weights(0,0,1) + weights(0,1,1) + weights(1,0,1) + weights(1,1,1) ! fractional position along x
                yf = weights(0,1,0) + weights(0,1,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along y
                zf = weights(1,0,0) + weights(1,0,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along z
                xfc=one-xf
                yfc=one-yf
                zfc=one-zf

                lim_x = dx(1) * xf * xfc
                xs = - prefactor * lim_x * (&
                            zfc * (&
                                    yfc * (phi(ks, js,   is+1) - phi(ks, js,   is))  &
                             +      yf *  (phi(ks, js+1, is+1) - phi(ks, js+1, is))) &
                     +      zf *  (&
                                    yfc * (phi(ks+1, js,   is+1) - phi(ks+1, js,   is))  &
                             +      yf *  (phi(ks+1, js+1, is+1) - phi(ks+1, js+1, is))))

                lim_x = lim_x * max_compression
                xs = max(-lim_x, min(xs, lim_x))

                lim_y = dx(2) * yf * yfc
                ys = - prefactor * lim_y * (&
                            zfc * (&
                                    xfc * (phi(ks, js+1, is  ) - phi(ks, js, is))    &
                             +      xf  * (phi(ks, js+1, is+1) - phi(ks, js, is+1))) &
                     +      zf *  (&
                                    xfc * (phi(ks+1, js+1, is) -   phi(ks+1, js, is))  &
                             +      xf  * (phi(ks+1, js+1, is+1) - phi(ks+1, js, is+1))))

                lim_y = lim_y * max_compression
                ys = max(-lim_y, min(ys, lim_y))

                lim_z = dx(3) * zf * zfc
                zs = - prefactor * lim_z * (&
                            xfc * (&
                                    yfc * (phi(ks+1, js, is  ) - phi(ks, js,   is))  &
                             +      yf  * (phi(ks+1, js+1, is) - phi(ks, js+1, is))) &
                     +      xf *  (&
                                    yfc * (phi(ks+1, js,   is+1) - phi(ks, js,   is+1))  &
                             +      yf  * (phi(ks+1, js+1, is+1) - phi(ks, js+1, is+1))))

                lim_z = lim_z * max_compression
                zs = max(-lim_z, min(zs, lim_z))

                parcels%position(1, n) = parcels%position(1, n) + xs
                parcels%position(2, n) = parcels%position(2, n) + ys
                parcels%position(3, n) = parcels%position(3, n) + zs
            enddo
            !$omp end do
            !$omp end parallel

            call parcel_communicate(parcels)

            call stop_timer(grad_corr_timer)

        end subroutine apply_gradient

end module parcel_correction
