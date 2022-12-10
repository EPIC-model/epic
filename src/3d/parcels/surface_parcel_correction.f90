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

    integer :: surf_lapl_corr_timer, &
               surf_grad_corr_timer

    private :: diverge

    public :: apply_surface_laplace,   &
              apply_surface_gradient,  &
              surf_lapl_corr_timer,    &
              surf_grad_corr_timer


    contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_surface_laplace
        call do_apply_laplace(lo_surf_parcels, n_lo_surf_parcels, 'lo')
        call do_apply_laplace(up_surf_parcels, n_up_surf_parcels, 'up')
    end subroutine apply_surface_laplace

    subroutine do_apply_laplace(s_parcels, n_par, which)
        type(surface_parcel_container_type), intent(inout) :: s_parcels
        integer,                             intent(in)    :: n_par
        character(2),                        intent(in)    :: which
        double precision                                   :: phi(0:ny-1,0:nx-1), &
                                                              ud(0:ny-1,0:nx-1),  &
                                                              vd(0:ny-1,0:nx-1)
        double precision                                   :: weights(ngp)
        integer                                            :: n, l, is(ngp), js(ngp), k

        call start_timer(surf_lapl_corr_timer)

        k = nz
        if (which == 'lo') then
            k = 0
        endif

        call do_area2grid(s_parcels, n_par, which)

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(k, :, :) * acelli - one

        call diverge(phi, ud, vd)

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud,vd) field:
        !$omp parallel default(shared)
        !$omp do private(n, l, is, js, weights)
        do n = 1, n_par
            call bilinear(s_parcels%position(:, n), is, js, weights)

            do l = 1, ngp
                s_parcels%position(1, n) = s_parcels%position(1, n)             &
                                         + weights(l) * ud(js(l), is(l))

                s_parcels%position(2, n) = s_parcels%position(2, n)             &
                                         + weights(l) * vd(js(l), is(l))
            enddo

            call apply_periodic_bc(s_parcels%position(:, n))
        enddo
        !$omp end do
        !$omp end parallel


        call stop_timer(surf_lapl_corr_timer)

    end subroutine do_apply_laplace

    subroutine apply_surface_gradient(prefactor, max_compression)
        double precision, intent(in) :: prefactor
        double precision, intent(in) :: max_compression
        call do_apply_gradient(lo_surf_parcels, n_lo_surf_parcels, 'lo', prefactor, max_compression)
        call do_apply_gradient(up_surf_parcels, n_up_surf_parcels, 'up', prefactor, max_compression)
    end subroutine apply_surface_gradient

    subroutine do_apply_gradient(s_parcels, n_par, which, prefactor, max_compression)
        type(surface_parcel_container_type), intent(inout) :: s_parcels
        integer,                             intent(in)    :: n_par
        character(2),                        intent(in)    :: which
        double precision,                    intent(in)    :: prefactor
        double precision,                    intent(in)    :: max_compression
        double precision                                   :: phi(0:ny-1,0:nx-1)
        double precision                                   :: weights(ngp)
        double precision                                   :: shift_x1, shift_x2
        double precision                                   :: x1_fpos, x2_fpos, lim_x1, lim_x2
        integer                                            :: n, is(ngp), js(ngp), k

        call start_timer(surf_grad_corr_timer)

        k = nz
        if (which == 'lo') then
            k = 0
        endif

        call do_area2grid(s_parcels, n_par, which)

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(k, :, :) * acelli - one

        !$omp parallel default(shared)
        !$omp do private(n, is, js, weights, x1_fpos, x2_fpos, shift_x1, shift_x2, lim_x1, lim_x2)
        do n = 1, n_par

            call bilinear(s_parcels%position(:, n), is, js, weights)

            x1_fpos=weights(2)+weights(4) ! fractional position along x1
            x2_fpos=weights(3)+weights(4) ! fractional position along x2

            shift_x1= - prefactor*dx(1)*x1_fpos*(one-x1_fpos)*(&
                        (one-x2_fpos)*(phi(js(2), is(2))-phi(js(1), is(1)))  &
                      +     (x2_fpos)*(phi(js(4), is(4))-phi(js(3), is(3))))

            lim_x1=max_compression*dx(1)*x1_fpos*(one-x1_fpos)

            shift_x1= max(-lim_x1,min(shift_x1,lim_x1))

            shift_x2= - prefactor*dx(2)*x2_fpos*(one-x2_fpos)*(&
                        (one-x1_fpos)*(phi(js(3), is(3))-phi(js(1), is(1))) &
                      +     (x1_fpos)*(phi(js(4), is(4))-phi(js(2), is(2))))

            lim_x2=max_compression*dx(2)*x2_fpos*(one-x2_fpos)

            shift_x2= max(-lim_x2,min(shift_x2,lim_x2))

            s_parcels%position(1, n) = s_parcels%position(1, n) + shift_x1
            s_parcels%position(2, n) = s_parcels%position(2, n) + shift_x2
        enddo
        !$omp end do
        !$omp end parallel

        call stop_timer(surf_grad_corr_timer)

    end subroutine do_apply_gradient

    ! Computes a divergent flow field (ud, vd) = grad(phi) where
    ! Lap(phi) = div (given).
    subroutine diverge(div, ud, vd)
        double precision, intent(inout) :: div(0:ny-1, 0:nx-1)
        double precision, intent(out)   :: ud(0:ny-1, 0:nx-1), vd(0:ny-1, 0:nx-1)
        double precision                :: ds(0:nx-1, 0:ny-1)
        double precision                :: us(0:nx-1, 0:ny-1), vs(0:nx-1, 0:ny-1)

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

    end subroutine

end module surface_parcel_correction
