!===================================================================
!  Uses a divergent flow to push parcels toward or away from grid
!  points so as to bring the area fraction at each grid point back
!  to unity.

!  Here tri-diagonal solves are used in the vertical direction.

!           Initialise data first using inipar.f90
!===================================================================

module parcel_diverge
    ! import FFT library:
    use stafft
    use deriv1d

    use constants
    use parameters, only : vcell

    use parcel_container

    implicit none

 !Gridded area fraction:
! c double precision:: volg(0:nz,nx)

    ! tri-diagonal arrays:
    double precision :: ap(nx),apb(nx)
    double precision :: etdv(0:nz,nx),htdv(0:nz,nx)
    double precision :: etd1(nz-1),htd1(nz-1)
    double precision :: etda(nz),htda(nz)

    ! wavenumbers::
    double precision :: hrkx(nx),rkx(nx)

    ! quantities needed in FFTs:
    double precision :: xtrig(2*nx)
    integer          :: xfactors(5)

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine init_diverge
            double precision:: a0(nx),a0b(nx),ksq(nx)
            double precision,parameter:: dzisq=dzi**2
            integer,parameter:: nwx=nx/2
            integer:: i,ix,kx,iz

            !----------------------------------------------------------
            ! set up FFTs:
            call initfft(nx,xfactors,xtrig)

            ! define x wavenumbers:
            call init_deriv(nx,ellx,hrkx)
            rkx(1) = zero
            do kx = 1, nwx-1
                rkx(kx+1)    = hrkx(2*kx)
                rkx(nx+1-kx) = hrkx(2*kx)
            enddo
            rkx(nwx+1) = hrkx(nx)

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
            htdv(:,1) = zero
            etdv(:,1) = zero
            do kx = 2,nx
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

    end subroutine init_diverge

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_diverge(volg)
        double precision, intent(in) :: volg(0:, -1:, :)
        double precision             :: phi(0:nz,nx),ud(0:nz,nx),wd(0:nz,nx)
        double precision             :: wbar(0:nz)
        double precision             :: weight(4)
        integer                      :: i, ngp, ij(2, 4)

        ! form divergence field * dt and store in phi temporarily:
        phi = volg - vcell

        !-----------------------------------------
        ! Forward x FFT:
        call forfft(nz+1,nx,phi,xtrig,xfactors)

        ! Compute the x-independent part of wd by integration:
        call vertint(phi(0,1),wbar)

        ! Invert Laplace's operator semi-spectrally with compact differences:
        call lapinv1(phi)

        ! Compute x derivative spectrally:
        call deriv(nz+1,nx,hrkx,phi,ud)

        ! Reverse x FFT to define x velocity component ud:
        call revfft(nz+1,nx,ud,xtrig,xfactors)

        ! Compute z derivative by compact differences:
        call diffz1(phi,wd)

        ! Add on the x-independent part of wd:
        wd(:,1) = wd(:,1) + wbar

        ! Reverse x FFT:
        call revfft(nz+1,nx,wd,xtrig,xfactors)

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud,wd) field:
        do i = 1, n_parcels
            call trilinear(parcels%position(n, :), ij, weight, ngp)

            do i = 1, ngp
                parcels%position(n, 1) = parcels%position(n, 1)             &
                                       + weight(i) * ud(ij(1, i)+1, ij(2, i))

                parcels%position(n, 2) = parcels%position(n, 2)             &
                                       + weight(i) * wd(ij(1, i)+1, ij(2, i))
            enddo

            call apply_periodic_bc(parcels%position(n, :))
        enddo
    end subroutine apply_diverge

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Inverts Laplace's operator on fs in semi-spectral space.
    ! Here dfs/dz = 0 on the z boundaries.
    ! Uses 4th-order compact differencing
    ! *** Overwrites fs ***
    subroutine lapinv1(fs)
        double precision, intent(inout) :: fs(0:nz,nx)
        double precision                :: rs(0:nz,nx)
        integer                         :: kx, iz

        fs(:,1) = zero
        do kx = 2, nx
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
        double precision, intent(in)  :: fs(0:nz,nx)
        double precision, intent(out) :: ds(0:nz,nx)
        integer                       :: ix, iz

        do ix = 1, nx
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

end module parcel_diverge
