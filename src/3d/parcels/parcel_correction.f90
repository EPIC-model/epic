!===================================================================
!       Module with divergent flow and gradient correction
!===================================================================
module parcel_correction
    use inversion_utils, only : init_fft
    use inversion_mod, only : diverge
    use parcel_interpl, only : trilinear, ngp, vol2grid
    use parcel_bc
    use omp_lib
    use constants
    use parameters, only : vcelli, nx, ny, nz, dx
    use parcel_container
    use timer, only : start_timer, stop_timer
    use fields, only : volg
    implicit none

    integer :: lapl_corr_timer, &
               grad_corr_timer

    public :: init_parcel_correction, &
              apply_laplace,          &
              apply_gradient,         &
              lapl_corr_timer,        &
              grad_corr_timer


    contains

        ! Initialise parcel correction module
        subroutine init_parcel_correction
            call init_fft
        end subroutine init_parcel_correction

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_laplace(l_reuse)
        logical, optional, intent(in) :: l_reuse
        double precision              :: phi(0:nz, 0:ny-1, 0:nx-1),   &
                                         ud(-1:nz+1, 0:ny-1, 0:nx-1), &
                                         vd(-1:nz+1, 0:ny-1, 0:nx-1), &
                                         wd(-1:nz+1, 0:ny-1, 0:nx-1)
        double precision              :: weights(ngp)
        integer                       :: n, l, is(ngp), js(ngp), ks(ngp)

        call start_timer(lapl_corr_timer)

        call vol2grid(l_reuse)

        ! form divergence field
        phi = volg(0:nz, :, :) * vcelli - one

        call diverge(phi, ud(0:nz, :, :), vd(0:nz, :, :), wd(0:nz, :, :))

        ! use symmetry to fill z grid lines outside domain:
        ud(-1, :, :)   =  ud(1, :, :)
        vd(-1, :, :)   =  vd(1, :, :)
        wd(-1, :, :)   = -wd(1, :, :)
        ud(nz+1, :, :) =  ud(nz-1, :, :)
        vd(nz+1, :, :) =  vd(nz-1, :, :)
        wd(nz+1, :, :) = -wd(nz-1, :, :)

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud, vd, wd) field:
        !$omp parallel default(shared)
        !$omp do private(n, l, is, js, ks, weights)
        do n = 1, n_parcels
            call trilinear(parcels%position(:, n), is, js, ks, weights)

            do l = 1, ngp
                parcels%position(1, n) = parcels%position(1, n)               &
                                       + weights(l) * ud(ks(l), js(l), is(l))

                parcels%position(2, n) = parcels%position(2, n)               &
                                       + weights(l) * vd(ks(l), js(l), is(l))

                parcels%position(3, n) = parcels%position(3, n)               &
                                       + weights(l) * wd(ks(l), js(l), is(l))
            enddo

            call apply_periodic_bc(parcels%position(:, n))
        enddo
        !$omp end do
        !$omp end parallel

        call stop_timer(lapl_corr_timer)

    end subroutine apply_laplace

    subroutine apply_gradient(prefactor, max_compression, l_reuse)
        double precision,  intent(in) :: prefactor
        double precision,  intent(in) :: max_compression
        logical, optional, intent(in) :: l_reuse
        double precision              :: phi(0:nz, 0:ny-1, 0:nx-1)
        double precision              :: weights(ngp)
        double precision              :: xs, ys, zs, xf, yf, zf, lim_x, lim_y, lim_z
        integer                       :: n, is(ngp), js(ngp), ks(ngp)

        call start_timer(grad_corr_timer)

        call vol2grid(l_reuse)

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(0:nz, :, :) * vcelli - one

        !$omp parallel default(shared)
        !$omp do private(n, is, js, ks, weights, xf, yf, zf, xs, ys, zs, lim_x, lim_y, lim_z)
        do n = 1, n_parcels

            call trilinear(parcels%position(:, n), is, js, ks, weights)

            xf = weights(2) + weights(4) + weights(6) + weights(8) ! fractional position along x
            yf = weights(3) + weights(4) + weights(7) + weights(8) ! fractional position along y
            zf = weights(5) + weights(6) + weights(7) + weights(8) ! fractional position along z

            ! l   ks(l) js(l) is(l)
            ! 1   k     j     i
            ! 2   k     j     i+1
            ! 3   k     j+1   i
            ! 4   k     j+1   i+1
            ! 5   k+1   j     i
            ! 6   k+1   j     i+1
            ! 7   k+1   j+1   i
            ! 8   k+1   j+1   i+1
            lim_x = dx(1) * xf * (one - xf)
            xs = - prefactor * lim_x * (&
                        (one - zf) * (&
                                (one - yf) * (phi(ks(2), js(2), is(2)) - phi(ks(1), js(1), is(1)))  &
                              +       (yf) * (phi(ks(4), js(4), is(4)) - phi(ks(3), js(3), is(3)))) &
                       +      (zf) * (&
                                (one - yf) * (phi(ks(6), js(6), is(6)) - phi(ks(5), js(5), is(5)))  &
                              +       (yf) * (phi(ks(8), js(8), is(8)) - phi(ks(7), js(7), is(7)))))

            lim_x = lim_x * max_compression
            xs = max(-lim_x, min(xs, lim_x))

            lim_y = dx(2) * yf * (one - yf)
            ys = - prefactor * lim_y * (&
                        (one - zf) * (&
                                (one - xf) * (phi(ks(3), js(3), is(3)) - phi(ks(1), js(1), is(1)))  &
                              +       (xf) * (phi(ks(4), js(4), is(4)) - phi(ks(2), js(2), is(2)))) &
                       +      (zf) * (&
                                (one - xf) * (phi(ks(7), js(7), is(7)) - phi(ks(5), js(5), is(5)))  &
                              +       (xf) * (phi(ks(8), js(8), is(8)) - phi(ks(6), js(6), is(6)))))

            lim_y = lim_y * max_compression
            ys = max(-lim_y, min(ys, lim_y))

            lim_z = dx(3) * zf * (one - zf)
            zs = - prefactor * lim_z * (&
                        (one - xf) * (&
                                (one - yf) * (phi(ks(5), js(5), is(5)) - phi(ks(1), js(1), is(1)))  &
                              +       (yf) * (phi(ks(7), js(7), is(7)) - phi(ks(3), js(3), is(3)))) &
                       +      (xf) * (&
                                (one - yf) * (phi(ks(6), js(6), is(6)) - phi(ks(2), js(2), is(2)))  &
                              +       (yf) * (phi(ks(8), js(8), is(8)) - phi(ks(4), js(4), is(4)))))

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
