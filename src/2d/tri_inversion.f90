! =============================================================================
!                       Inversion of 2D vorticity
!               (the component out of the plane of the flow)
!
!          Here tri-diagonal solves are used in the vertical direction.
! =============================================================================

module tri_inversion
    use stafft
    use deriv1d
    use constants
    use parameters, only : nx, nz, dx, dxi, extent, lower
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,  &
               vtend_timer

    private
        ! Tri-diagonal arrays:
        double precision, allocatable :: ap(:)
        double precision, allocatable :: etdh(:, :), htdh(:, :)
        double precision, allocatable :: etd0(:), htd0(:)

        ! Wavenumbers::
        double precision, allocatable :: hrkx(:), rkx(:)

        ! Quantities needed in FFTs:
        double precision, allocatable :: xtrig(:)
        integer                       :: xfactors(5)

        logical :: l_filled = .false.

    public :: init_inversion,       &
              vor2vel,              &
              vorticity_tendency,   &
              vor2vel_timer,        &
              vtend_timer

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            double precision :: a0(0:nx-1), ksq(0:nx-1)
            double precision :: dzisq
            integer          :: nwx, kx, iz

            dzisq = dxi(2) ** 2
            nwx = nx / 2

            if (.not. allocated(ap)) then
                allocate(ap(0:nx-1))
                allocate(etdh(nz-1,0:nx-1))
                allocate(htdh(nz-1,0:nx-1))
                allocate(etd0(0:nz))
                allocate(htd0(0:nz))
                allocate(hrkx(0:nx-1))
                allocate(rkx(0:nx-1))
                allocate(xtrig(2*nx))
            endif

            !----------------------------------------------------------
            ! Set up FFTs:
            call initfft(nx, xfactors, xtrig)

            ! Define x wavenumbers:
            call init_deriv(nx, extent(1), hrkx)
            rkx(0) = zero
            do kx = 1, nwx-1
                rkx(kx)    = hrkx(2*kx-1)
                rkx(nx-kx) = hrkx(2*kx-1)
            enddo
            rkx(nwx) = hrkx(nx-1)

            ! Squared wavenumber array (used in tridiagonal solve):
            ksq = rkx ** 2

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0 = -two * dzisq - f56 * ksq
            ap = dzisq - f112*ksq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for inversion of Poisson's equation (Dirichlet BC):
            do kx = 1,nx-1
                htdh(1,kx) = one / a0(kx)
                etdh(1,kx) = -ap(kx) * htdh(1,kx)
                do iz = 2,nz-2
                    htdh(iz,kx) = one / (a0(kx) + ap(kx) * etdh(iz-1,kx))
                    etdh(iz,kx) = -ap(kx) * htdh(iz,kx)
                enddo
                htdh(nz-1,kx) = one / (a0(kx) + ap(kx) * etdh(nz-2,kx))
            enddo

            ! Tridiagonal arrays for the compact difference calculation of d/dz
            ! for fields f for which f = 0 at the boundaries:
            htd0(0) = one / f23
            etd0(0) = -f13 * htd0(0)
            do iz = 1,nz-1
                htd0(iz) = one / (f23 + f16 * etd0(iz-1))
                etd0(iz) = -f16 * htd0(iz)
            enddo
            htd0(nz) = one / (f23 + f13 * etd0(nz-1))
        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Inverts the vorticity "vortg" to obtain the gridded velocity field
        ! u = velog(:, :, 1) = -dpsig/dz and w = velog(:, :, 2) = dpsig/dx
        ! and computes the velocity gradient "velgradg".
        subroutine vor2vel(vortg, velog, velgradg)
            double precision, intent(in)  :: vortg(-1:nz+1, 0:nx-1)
            double precision, intent(out) :: velog(-1:nz+1, 0:nx-1, 2)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:nx-1, 4)
            double precision              :: xx, zz, sinx, sinz, cosx, cosz
            integer                       :: iz, ix


            if (l_filled) then
                return
            endif

            l_filled = .true.

            call start_timer(vor2vel_timer)

            do ix = 0, nx - 1
                do iz = -1, nz+1
                    xx = pi * (lower(1) + dx(1) * dble(ix))
                    zz = pi * (lower(2) + dx(2) * dble(iz))

                    sinx = dsin(xx)
                    cosx = dcos(xx)

                    sinz = dsin(zz)
                    cosz = dcos(zz)

                    ! horizontal velocity, u
                    velog(iz, ix, 1) = sinx * cosz

                    ! vertical velocity, w
                    velog(iz, ix, 2) = -cosx * sinz

                    ! du/dx
                    velgradg(iz, ix, 1) = pi * cosx * cosz

                    ! dw/dx
                    velgradg(iz, ix, 3) = pi * sinx * sinz

                    ! du/dz = - dw/dx
                    velgradg(iz, ix, 2) = - velgradg(iz, ix, 3)

                    ! dw/dz = - du/dx
                    velgradg(iz, ix, 4) = - velgradg(iz, ix, 1)

                enddo
            enddo

            call stop_timer(vor2vel_timer)

        end subroutine vor2vel


        subroutine vorticity_tendency(tbuoyg, vtend)
            double precision, intent(in)  :: tbuoyg(-1:nz+1, 0:nx-1)
            double precision, intent(out) :: vtend(-1:nz+1, 0:nx-1)
            double precision              :: psig(0:nz, 0:nx-1)

            call start_timer(vtend_timer)

            psig = tbuoyg(0:nz, 0:nx-1)

            ! Forward x FFT:
            call forfft(nz+1, nx, psig, xtrig, xfactors)

            ! Compute x derivative spectrally of psig:
            call deriv(nz+1, nx, hrkx, psig, vtend(0:nz, :))

            ! Reverse x FFT
            call revfft(nz+1, nx, vtend(0:nz, :), xtrig, xfactors)

            ! Fill z grid lines outside domain:
            vtend(-1,   :) = two * vtend(0,  :) - vtend(1,    :)
            vtend(nz+1, :) = two * vtend(nz, :) - vtend(nz-1, :)

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here fs = 0 on the z boundaries.
        !Uses 4th-order compact differencing
        !*** Overwrites fs ***
        subroutine lapinv0(fs)
            double precision, intent(inout) :: fs(0:nz,0:nx-1)
            double precision                :: rs(nz-1,0:nx-1)
            integer                         :: kx,iz

            fs(:,0) = zero !x-independent mode cannot be inverted this way
            do kx = 1,nx-1
                do iz = 1,nz-1
                    rs(iz,kx) = f112 * (fs(iz-1,kx) + fs(iz+1,kx)) + f56 * fs(iz,kx)
                enddo

                fs(0,kx) = zero
                fs(1,kx) = rs(1,kx) * htdh(1,kx)
                do iz =2,nz-1
                    fs(iz,kx) = (rs(iz,kx) - ap(kx) * fs(iz-1,kx)) * htdh(iz,kx)
                enddo
                fs(nz,kx) = zero

                do iz = nz-2,1,-1
                    fs(iz,kx) = etdh(iz,kx) * fs(iz+1,kx) + fs(iz,kx)
                enddo
            enddo
        end subroutine lapinv0

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f which has f = 0 at the boundaries
        !using 4th-order compact differencing.  Here fs = f, ds = df/dz,
        !lapfsbot = Laplace(f) at z_min and lapfstop = Laplace(f) at z_max,
        !both given.  *** All quantities must be in semi-spectral space ***
        subroutine diffz0(fs, ds, lapfsbot, lapfstop)
            double precision, intent(in)  :: fs(0:nz,0:nx-1)
            double precision, intent(out) :: ds(0:nz,0:nx-1)
            double precision, intent(in)  :: lapfsbot(0:nx-1), lapfstop(0:nx-1)
            double precision              :: hdzi, dz6
            integer                       :: kx, iz

            dz6  = f16 * dx(2)
            hdzi = f12 * dxi(2)

            do kx = 0,nx-1
                ds(0,kx) = fs(1,kx) * dxi(2) - dz6 * lapfsbot(kx)
                ds(1,kx) = fs(2,kx) * hdzi
                do iz = 2,nz-2
                    ds(iz,kx) = (fs(iz+1,kx) - fs(iz-1,kx)) * hdzi
                enddo
                ds(nz-1,kx) = -fs(nz-2,kx) * hdzi
                ds(nz,kx) = dz6 * lapfstop(kx) - fs(nz-1,kx) * dxi(2)

                ds(0,kx) = ds(0,kx) * htd0(0)
                do iz = 1,nz-1
                    ds(iz,kx) = (ds(iz,kx) - f16 * ds(iz-1,kx)) * htd0(iz)
                enddo
                ds(nz,kx) = (ds(nz,kx) - f13 * ds(nz-1,kx)) * htd0(nz)

                do iz = nz-1,0,-1
                    ds(iz,kx) = etd0(iz) * ds(iz+1,kx) + ds(iz,kx)
                enddo
            enddo
        end subroutine diffz0
end module tri_inversion
