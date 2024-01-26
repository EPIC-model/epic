!===================================================================
!  Uses a divergent flow to push parcels toward or away from grid
!  points so as to bring the area fraction at each grid point back
!  to unity.

!  Here tri-diagonal solves are used in the vertical direction.

!           Initialise data first using inipar.f90
!===================================================================

module surface_parcel_correction
    use sta2dfft, only : ptospc, spctop, xderiv, yderiv
    use inversion_utils, only : lapinv, hrkx, hrky, xfactors, yfactors, xtrig, ytrig
    use surface_parcel_interpl, only : bilinear, ngp, do_area2grid
    use parcel_bc
    use omp_lib
    use constants
    use parameters, only : acelli, nx, ny, dx, dxi, nz
    use surface_parcel_container
    use timer, only : start_timer, stop_timer
    use fields, only : volg

    implicit none

    private

!     integer :: surf_lapl_corr_timer, &
!                surf_grad_corr_timer

    public :: apply_surf_laplace,   &
              apply_surf_gradient!,  &
!               surf_lapl_corr_timer,    &
!               surf_grad_corr_timer


    contains


    subroutine apply_surf_laplace

        call area2grid(l_halo_swap=.true.)

        call m_apply_laplace(0,  bot_parcels)
        call m_apply_laplace(nz, top_parcels)
    end subroutine apply_surf_laplace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_apply_laplace(iz, surf_parcels)
        integer,               intent(in)    :: iz
        type(ellipse_pc_type), intent(inout) :: surf_parcels
        double precision                     :: phi( box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision                     :: grad(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 2)
        double precision                     :: weights(0:1, 0:1, 0:1)
        integer                              :: n, l, is, js

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(iz, :, :) * acelli - one

        call m_diverge(phi, grad)

        !------------------------------------------------------------------
        ! Increment parcel positions using (ud, vd) field:
        !$omp parallel default(shared)
        !$omp do private(n, l, is, js, weights)
        do n = 1, n_par
            call bilinear(surf_parcels%position(:, n), is, js, weights)

            do l = 1, 2
                surf_parcels%position(1, n) = surf_parcels%position(1, n)               &
                                            + sum(weights * grad(js:js+1, is:is+1, l))
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        !FIXME
        call parcel_communicate

    end subroutine m_apply_laplace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_surf_gradient(prefactor, max_compression)
        double precision, intent(in) :: prefactor
        double precision, intent(in) :: max_compression

        call area2grid(l_halo_swap=.true.)

        call m_apply_gradient(0,  bot_parcels, prefactor, max_compression)
        call m_apply_gradient(nz, top_parcels, prefactor, max_compression)

    end subroutine apply_surf_gradient

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_apply_gradient(iz, surf_parcels, prefactor, max_compression)
        integer,               intent(in)    :: iz
        type(ellipse_pc_type), intent(inout) :: surf_parcels
        double precision,      intent(in)    :: prefactor
        double precision,      intent(in)    :: max_compression
        double precision                     :: phi(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision                     :: weights(0:1, 0:1)
        double precision                     :: xs, ys, xf, yf, xfc, yfc, lim_x, lim_y
        integer                              :: n, is, js

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(iz, :, :) * acelli - one

        !$omp parallel default(shared)
        !$omp do private(n, is, js, weights, x1_fpos, x2_fpos, shift_x1, shift_x2, lim_x1, lim_x2)
        do n = 1, surf_parcels%local_num

            call bilinear(surf_parcels%position(:, n), is, js, weights)

            xf = weights(0, 1) + weights(1, 1) ! fractional position along x
            yf = weights(1, 0) + weights(1, 1) ! fractional position along y
            xfc = one - xf
            yfc = one - yf


            lim_x = dx(1) * xf * xfc
            xs = - prefactor * lim_x * (                             &
                            yfc * (phi(js,   is+1) - phi(js,   is))  &
                          + yf  * (phi(js+1, is+1) - phi(js+1, is)))

            lim_x = lim_x * max_compression
            xs = max(-lim_x, min(xs, lim_x))

            lim_y = dx(2) * yf * yfc
            ys = - prefactor * lim_y * (                            &
                            xfc * (phi(js+1, is  ) - phi(js, is))   &
                          + xf  * (phi(js+1, is+1) - phi(js, is+1)))

            lim_y = lim_y * max_compression
            ys = max(-lim_y, min(ys, lim_y))

            surf_parcels%position(1, n) = surf_parcels%position(1, n) + xs
            surf_parcels%position(2, n) = surf_parcels%position(2, n) + ys
        enddo
        !$omp end do
        !$omp end parallel

        !FIXME
        call parcel_communicate

    end subroutine m_apply_gradient

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !FIXME
    ! Computes a divergent flow field (ud, vd) = grad(phi) where
    ! Lap(phi) = div (given).
    subroutine m_diverge(div, grad)
        double precision, intent(inout) :: div(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision, intent(out)   :: grad(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 2)
        double precision                :: ds(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                :: us(box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                                           vs(box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        !------------------------------------------------------------------
        ! Convert phi to spectral space as ds:
        call ptospc(nx, ny, div, ds, xfactors, yfactors, xtrig, ytrig)

        ! Invert Laplace's operator spectrally
        call lapinv(ds)

        ! Compute x derivative spectrally:
        call xderiv(nx, ny, hrkx, ds, us)

        ! Reverse FFT to define x velocity component ud:
        call spctop(nx, ny, us, ud, xfactors, yfactors, xtrig, ytrig)

        ! Compute y derivative spectrally:
        call yderiv(nx, ny, hrky, ds, vs)

        ! Reverse FFT to define y velocity component vd:
        call spctop(nx, ny, vs, vd, xfactors, yfactors, xtrig, ytrig)

    end subroutine m_diverge

end module surface_parcel_correction
