!===================================================================
!  Uses a divergent flow to push parcels toward or away from grid
!  points so as to bring the area fraction at each grid point back
!  to unity.

!  Here tri-diagonal solves are used in the vertical direction.

!           Initialise data first using inipar.f90
!===================================================================

module parcel_correction
    ! import FFT library:
    use stafft
    use deriv1d

    use parcel_interpl, only : bilinear, ngp, vol2grid
    use parcel_bc
    use omp_lib

    use constants
    use parameters, only : vcelli, nx, nz, dx, dxi, lcelli

    use parcel_container

    use timer, only : start_timer, stop_timer

    use fields, only : volg

    implicit none

    integer :: lapl_corr_timer, &
               grad_corr_timer

    private
        ! tri-diagonal arrays:
        double precision, allocatable :: ap(:), apb(:)
        double precision, allocatable :: etdv(:, :), htdv(:, :)
        double precision, allocatable :: etd1(:), htd1(:)
        double precision, allocatable :: etda(:), htda(:)

        ! wavenumbers::
        double precision, allocatable :: hrkx(:), rkx(:)

        ! quantities needed in FFTs:
        double precision, allocatable :: xtrig(:)
        integer                       :: xfactors(5)

    public :: init_parcel_correction, &
              apply_laplace,          &
              apply_gradient,         &
              lapl_corr_timer,        &
              grad_corr_timer


    contains

        ! Initialise parcel correction module. It allocates memory and does
        ! some precomputations.
        subroutine init_parcel_correction
            double precision :: a0(0:nx-1), a0b(0:nx-1), ksq(0:nx-1)
            double precision :: dzisq
            integer          :: nwx, kx, iz

            dzisq = dxi(2) ** 2
            nwx = nx / 2

            if (.not. allocated(ap)) then
                allocate(ap(0:nx-1))
                allocate(apb(0:nx-1))
                allocate(etdv(0:nz,0:nx-1))
                allocate(htdv(0:nz,0:nx-1))
                allocate(etd1(nz-1))
                allocate(htd1(nz-1))
                allocate(etda(nz))
                allocate(htda(nz))
                allocate(hrkx(0:nx-1))
                allocate(rkx(0:nx-1))
                allocate(xtrig(2*nx))
            endif

            !----------------------------------------------------------
            ! set up FFTs:
            call initfft(nx, xfactors, xtrig)

            ! define x wavenumbers:
            call init_deriv(nx, extent(1), hrkx)
            rkx(0) = zero
            do kx = 1, nwx-1
                rkx(kx)    = hrkx(2*kx-1)
                rkx(nx-kx) = hrkx(2*kx-1)
            enddo
            rkx(nwx) = hrkx(nx-1)

            ! squared wavenumber array (used in tridiagonal solve):
            ksq = rkx**2

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0  = - two * dzisq - f56 * ksq
            a0b =       - dzisq - f13 * ksq
            ap  =         dzisq - f112 * ksq
            apb =         dzisq - f16 * ksq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for inversion of Poisson's equation:
            htdv(:,0) = zero
            etdv(:,0) = zero
            do kx = 1,nx-1
                htdv(0,kx) = one / a0b(kx)
                etdv(0,kx) = -apb(kx) * htdv(0,kx)
                do iz = 1,nz-1
                    htdv(iz,kx) = one / (a0(kx) + ap(kx) * etdv(iz-1,kx))
                    etdv(iz,kx) = -ap(kx) * htdv(iz,kx)
                enddo
                htdv(nz,kx) = one / (a0b(kx) + apb(kx) * etdv(nz-1,kx))
            enddo

            ! Tridiagonal arrays for the compact difference calculation of d/dz
            ! for fields f for which df/dz = 0 at the boundaries:
            htd1(1) = one / f23
            etd1(1) = -f16 * htd1(1)
            do iz = 2,nz-2
                htd1(iz) = one / (f23 + f16 * etd1(iz-1))
                etd1(iz) = -f16 * htd1(iz)
            enddo
            htd1(nz-1) = one / (f23 + f16 * etd1(nz-2))

            ! Tridiagonal arrays used for integrating in z (see vertint):
            htda(1) = one / f76
            etda(1) = -f16 * htda(1)
            do iz = 2,nz-1
                htda(iz) = one / (one + f16 * etda(iz-1))
                etda(iz) = -f16 * htda(iz)
            enddo
            htda(nz) = one / (f76 + f16 * etda(nz-2))

    end subroutine init_parcel_correction

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_laplace
        double precision :: phi(0:nz,0:nx-1), ud(-1:nz+1,0:nx-1), wd(-1:nz+1,0:nx-1)
        double precision :: wbar(0:nz)
        double precision :: weights(ngp)
        integer          :: n, l, is(ngp), js(ngp)

        call start_timer(lapl_corr_timer)

        call vol2grid

        ! form divergence field * dt and store in phi temporarily:
        phi(1:nz-1, :) = volg(1:nz-1, :) * vcelli - one
        phi(0,      :) = volg(0,  :) * lcelli - one
        phi(nz,     :) = volg(nz, :) * lcelli - one

        !-----------------------------------------
        ! Forward x FFT:
        call forfft(nz+1,nx,phi,xtrig,xfactors)

        ! Compute the x-independent part of wd by integration:
        call vertint(phi(0,0),wbar)

        ! Invert Laplace's operator semi-spectrally with compact differences:
        call lapinv1(phi)

        ! Compute x derivative spectrally:
        call deriv(nz+1,nx,hrkx,phi,ud(0:nz, :))

        ! Reverse x FFT to define x velocity component ud:
        call revfft(nz+1,nx,ud(0:nz, :),xtrig,xfactors)

        ! Compute z derivative by compact differences:
        call diffz1(phi,wd(0:nz, :))

        ! Add on the x-independent part of wd:
        wd(0:nz,0) = wd(0:nz,0) + wbar

        ! Reverse x FFT:
        call revfft(nz+1,nx,wd(0:nz,:),xtrig,xfactors)

        ! Use symmetry to fill z grid lines outside domain:
        ud(-1,:)   =  ud(1,:)
        wd(-1,:)   = -wd(1,:)
        ud(nz+1,:) =  ud(nz-1,:)
        wd(nz+1,:) = -wd(nz-1,:)

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud,wd) field:
        !$omp parallel default(shared)
        !$omp do private(n, l, is, js, weights)
        do n = 1, n_parcels
            call bilinear(parcels%position(:, n), is, js, weights)

            do l = 1, ngp
                parcels%position(1, n) = parcels%position(1, n)             &
                                       + weights(l) * ud(js(l), is(l))

                parcels%position(2, n) = parcels%position(2, n)             &
                                       + weights(l) * wd(js(l), is(l))
            enddo

            call apply_periodic_bc(parcels%position(:, n))
        enddo
        !$omp end do
        !$omp end parallel


        call stop_timer(lapl_corr_timer)

    end subroutine apply_laplace

    subroutine apply_gradient(prefactor, max_compression)
        double precision, intent(in) :: prefactor
        double precision, intent(in) :: max_compression
        double precision             :: phi(0:nz,0:nx-1)
        double precision             :: weights(ngp)
        double precision             :: shift_x1, shift_x2, x1_fpos, x2_fpos, lim_x1, lim_x2
        integer                      :: n, is(ngp), js(ngp)

        call start_timer(grad_corr_timer)

        call vol2grid

        ! form divergence field * dt and store in phi temporarily:
        phi(1:nz-1, :) = volg(1:nz-1, :) * vcelli - one
        phi(0,      :) = volg(0,  :) * lcelli - one
        phi(nz,     :) = volg(nz, :) * lcelli - one

        !$omp parallel default(shared)
        !$omp do private(n, is, js, weights, x1_fpos, x2_fpos, shift_x1, shift_x2, lim_x1, lim_x2)
        do n = 1, n_parcels

            call bilinear(parcels%position(:, n), is, js, weights)

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

            parcels%position(1, n) = parcels%position(1, n) + shift_x1
            parcels%position(2, n) = parcels%position(2, n) + shift_x2
        enddo
        !$omp end do
        !$omp end parallel

        call stop_timer(grad_corr_timer)

    end subroutine apply_gradient
    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Inverts Laplace's operator on fs in semi-spectral space.
    ! Here dfs/dz = 0 on the z boundaries.
    ! Uses 4th-order compact differencing
    ! *** Overwrites fs ***
    subroutine lapinv1(fs)
        double precision, intent(inout) :: fs(0:nz,0:nx-1)
        double precision                :: rs(0:nz,0:nx-1)
        integer                         :: kx, iz

        fs(:,0) = zero
        do kx = 1, nx-1
            rs(0,kx) = f13 * fs(0,kx) + f16 * fs(1,kx)
            do iz = 1,nz-1
                rs(iz,kx) = f112 * (fs(iz-1,kx) + fs(iz+1,kx)) + f56 * fs(iz,kx)
            enddo
            rs(nz,kx) = f13 * fs(nz,kx) + f16 * fs(nz-1,kx)

            fs(0,kx) = rs(0,kx) * htdv(0,kx)
            do iz = 1, nz-1
                fs(iz,kx) = (rs(iz,kx) - ap(kx) * fs(iz-1,kx)) * htdv(iz,kx)
            enddo
            fs(nz,kx) = (rs(nz,kx) - apb(kx) * fs(nz-1,kx)) * htdv(nz,kx)

            do iz = nz-1, 0, -1
                fs(iz,kx) = etdv(iz,kx) * fs(iz+1,kx) + fs(iz,kx)
            enddo
        enddo
    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates df/dz for a field f which has df/dz = 0 at the boundaries
    ! using 4th-order compact differencing.  Here fs = f and ds = df/dz.
    subroutine diffz1(fs,ds)
        double precision, intent(in)  :: fs(0:nz,0:nx-1)
        double precision, intent(out) :: ds(0:nz,0:nx-1)
        integer                       :: ix, iz
        double precision              :: hdzi

        hdzi = f12 * dxi(2)

        do ix = 0, nx-1
            do iz = 1, nz-1
                ds(iz,ix) = (fs(iz+1,ix) - fs(iz-1,ix)) * hdzi
            enddo

            ds(0,ix) = zero
            ds(1,ix) = ds(1,ix) * htd1(1)
            do iz = 2, nz-1
                ds(iz,ix) = (ds(iz,ix) - f16 * ds(iz-1,ix)) * htd1(iz)
            enddo
            ds(nz,ix) = zero

            do iz = nz-2, 1, -1
                ds(iz,ix) = etd1(iz) * ds(iz+1,ix) + ds(iz,ix)
            enddo
        enddo
    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Finds f by integrating df/dz = d, ensuring f = 0 at the boundaries
    ! using 4th-order compact differencing.  Here ds = df/dz and fs = f.
    subroutine vertint(ds, fs)
        double precision, intent(in)  :: ds(0:nz)
        double precision, intent(out) :: fs(0:nz)
        double precision              :: es(nz), esum
        integer                       :: iz
        double precision              :: dz2
        double precision              :: hdzi
        double precision              :: dz24

        dz2  = f12 * dx(2)
        hdzi = f12 * dxi(2)
        dz24 = f124 * dx(2)

        !-------------------------------------------
        ! First interpolate ds to a half grid as es:
        do iz = 1, nz
            es(iz) = f23 * (ds(iz-1) + ds(iz))
        enddo

        es(1) = es(1) * htda(1)
        do iz = 2, nz
            es(iz) = (es(iz) - f16 * es(iz-1)) * htda(iz)
        enddo

        do iz = nz-1, 1, -1
            es(iz) = etda(iz) * es(iz+1) + es(iz)
        enddo

        !-------------------------------------------
        ! Next adjust es to ensure f(nz) = 0:
        esum = (f1112 * (es(1) + es(nz)) + sum(es(2:nz-1))) / (dble(nz) - f16)
        es = es - esum

        !Integrate:
        fs(0) = zero
        fs(1) = dz24 * (23.d0 * es(1) + es(2))
        do iz = 2, nz-1
            fs(iz) = fs(iz-1) + dz24 * (es(iz-1) + 22.d0 * es(iz) + es(iz+1))
        enddo
        fs(nz) = zero

    end subroutine

end module parcel_correction
