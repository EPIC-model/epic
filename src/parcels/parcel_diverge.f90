! =============================================================================
!       Uses a divergent flow to push parcels toward or away from grid
!       points so as to bring the area fraction at each grid point back
!       to unity.
! =============================================================================
module parcel_diverge

    ! import FFT library:
    use stafft
    use deriv1d

    use parameters, only : nx, nz, extent, vcell
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use interpl_methods, only : trilinear

    ! import parameters and constants:
    use constants, only : one, zero

    implicit none

    ! Wavenumbers and inverse Laplacian:
    double precision, allocatable :: hrkx(:), rkx(:), rkz(:)
    double precision, allocatable :: laplinv(:, :)

    ! Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ztrig(:)
    integer          :: xfactors(5), zfactors(5)

    contains

    ! Initialise FFT
    subroutine init_diverge
        integer :: i, kx, kz

        allocate(hrkx(nx))
        allocate(rkx(nx))
        allocate(rkz(0:nz+1))
        allocate(xtrig(2 * nx))
        allocate(ztrig(2 * (nz+2)))
        allocate(laplinv(-1:nz+1, nx))

        ! Set up FFTs:
        call initfft(nx, xfactors, xtrig)
        call initfft(nz+2, zfactors, ztrig)

        ! Define x wavenumbers:
        call init_deriv(nx, extent(1), hrkx)
        rkx(1) = zero
        do kx = 1, nx/2 - 1
            rkx(kx+1)    = hrkx(2*kx)
            rkx(nx+1-kx) = hrkx(2*kx)
        enddo
        rkx(nx/2 + 1) = hrkx(nx)

        ! Define z wavenumbers:
        call init_deriv(nz+2, extent(2), rkz)

        ! Define spectral inverse Laplacian for inverting Poisson's equation:
        do kx = 1, nx
            do kz = 0, nz+1
                laplinv(kz,kx) = -one / (rkx(kx)**2 + rkz(kz)**2)
            enddo
        enddo
        ! kz = 0:
        do kx=2,nx
            laplinv(-1, kx) = -one / rkx(kx)**2
        enddo
        ! The zero wavenumber mode has no significance:
        laplinv(-1, 1) = zero
        ! For z derivatives of a cosine in z function:
        rkz(nz+1) = zero

    end subroutine init_diverge


    subroutine apply_diverge(volg)
        double precision, intent(in) :: volg(0:, -1:, :)
        double precision             :: ud(nx, -1:nz+1),  wd(nx, -1:nz+1),  wka(nx, 0:nz+1)
        double precision             :: phi(-1:nz+1, nx), uds(-1:nz+1, nx), wds(0:nz+1, nx)
        integer                      :: i, n, ngp, ij(2, 4)
        integer                      :: ix, iz, kx, kz
        double precision             :: weight(4)
        double precision             :: pos(2)

        ! Form divergence field * dt and store in ud temporarily:
        ! (normalize volg by cell volume)
        ud = volg(:, :, 1) - vcell

        !-----------------------------------------
        ! Forward z cosine FFT:
        call dct(nx, nz+2, ud, ztrig, zfactors)
        ! Transpose array:
        do ix = 1, nx
            do kz = -1, nz+1
                phi(kz, ix) = ud(ix, kz)
            enddo
        enddo
        ! Forward x FFT:
        call forfft(nz+3, nx, phi, xtrig, xfactors)

        ! Invert Laplace's operator spectrally:
        phi = laplinv * phi
        ! phi = (spectral) velocity potential * dt

        !-----------------------------------------
        ! Compute x derivative spectrally:
        call deriv(nz+3, nx, hrkx, phi, uds)

        ! Reverse x FFT:
        call revfft(nz+3, nx, uds, xtrig, xfactors)
        ! Transpose array:
        do kz = -1, nz+1
            do ix = 1, nx
                ud(ix, kz) = uds(kz, ix)
            enddo
        enddo
        ! Reverse z cosine FFT:
        call dct(nx, nz+2, ud, ztrig, zfactors)

        !-----------------------------------------
        ! Compute z derivative spectrally:
        do kx = 1, nx
            do kz = 0, nz+1
                wds(kz, kx) = -rkz(kz) * phi(kz, kx)
            enddo
        enddo
        ! This makes wds a sine series in z

        ! Reverse x FFT:
        call revfft(nz+2, nx, wds, xtrig, xfactors)
        ! Transpose array:
        do kz = 0, nz+1
            do ix = 1, nx
                wka(ix, kz) = wds(kz, ix)
            enddo
        enddo
        ! Reverse z sine FFT:
        call dst(nx, nz+2, wka, ztrig, zfactors)

        ! Copy into wd with zero edge values:
        wd(:, -1) = zero
        do iz = 0, nz
            wd(:, iz) = wka(:, iz)
        enddo
        wd(:, nz+1) = zero

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud,wd) field:
        do n = 1, n_parcels

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

end module parcel_diverge
