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

    ! import parameters and constants:
    use constants, only : one, zero

    implicit none

!  !Total number of parcels:
! c integer:: n

!  !Parcel positions (x,z) and area fractions:
! c double precision:: x(nm),z(nm),v(nm)

!  !Rates of change (dx/dt,dz/dt):
! c double precision:: dxdt(nm),dzdt(nm)

!     !Gridded area fraction:
!     double precision:: volg(nx,0:nz)

    ! Wavenumbers and inverse Laplacian:
    double precision, allocatable :: hrkx(:), rkx(:), rkz(:)
    double precision, allocatable :: laplinv(:, :)

    ! Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ztrig(:)
    integer          :: xfactors(5), zfactors(5)

!  !Next horizontal grid point used in bilinear interpolation:
! c integer:: ixp(nx)

!  !Time step (does not matter):
! c double precision,parameter:: dt=0.01d0

! c double precision:: t
! c integer:: istep,nsteps

! !---------------------------------------------------------------------
! ! Initialise (read data, define fixed arrays):
! call initialise
!
! write(*,*) ' Enter the number of time steps used to relax the'
! write(*,*) ' gridded area fraction back to unity:'
! read(*,*) nsteps

! ! Evolve parcels:
! do istep=1,nsteps
!   call evolve(istep)
! enddo

! ! Finalise:
! call par2grid
!
! ! Write final data
! write(80,rec=nsteps+1) real(nsteps+1),real(volg-one)
!
! close(40)
! close(80)
!
! write(*,*)
! write(*,*) ' The rms area fraction error is in volf.asc and the'
! write(*,*) ' gridded area fraction error is in dv.r4; view with'
! write(*,*)
! write(*,'(a,i3,1x,i3,a)') '  dataview -ndim ',nx,nz+1,' dv.r4 -mod &'
! write(*,*)

    contains

    ! Initialise FFT
    subroutine init_diverge
        integer :: i, kx, kz

        allocate(hrkx(nx))
        allocate(rkx(nx))
        allocate(rkz(nz))
        allocate(xtrig(2 * nx))
        allocate(ztrig(2 * nz))
        allocate(laplinv(0:nz, nx))

        ! Set up FFTs:
        call initfft(nx, xfactors, xtrig)
        call initfft(nz, zfactors, ztrig)

        ! Define x wavenumbers:
        call init_deriv(nx, extent(1), hrkx)
        rkx(1) = zero
        do kx = 1, nx/2 - 1
            rkx(kx+1)    = hrkx(2*kx)
            rkx(nx+1-kx) = hrkx(2*kx)
        enddo
        rkx(nx/2 + 1) = hrkx(nx)

        ! Define z wavenumbers:
        call init_deriv(nz, extent(2), rkz)

        ! Define spectral inverse Laplacian for inverting Poisson's equation:
        do kx = 1, nx
            do kz = 1, nz
                laplinv(kz,kx) = -one / (rkx(kx)**2 + rkz(kz)**2)
            enddo
        enddo
        ! kz = 0:
        do kx=2,nx
            laplinv(0, kx) = -one / rkx(kx)**2
        enddo
        ! The zero wavenumber mode has no significance:
        laplinv(0, 1) = zero
        ! For z derivatives of a cosine in z function:
        rkz(nz) = zero

        ! !---------------------------------------------------------------
        ! ! File to contain rms area fraction - 1 & max abs value vs time:
        ! open(40,file='volf.asc',status='replace')
        !
        ! ! File to contain the gridded area fraction - 1 vs time:
        ! open(80,file='dv.r4',form='unformatted',access='direct', &
        !                    status='replace',recl=4*(nx*(nz+1)+1))
    end subroutine init_diverge


    subroutine apply_divergent_flow(volg)
        double precision, intent(in) :: volg(:, :)
        double precision             :: ud(nx, 0:nz),  wd(nx, 0:nz),  wka(nx, nz)
        double precision             :: phi(0:nz, nx), uds(0:nz, nx), wds(nz, nx)
        integer                      :: i, n, ngp, ij(2, 4)
        integer                      :: ix, iz, kx, kz
        double precision             :: weight(4)

! ! ! ! !-----------------------------------------------------------------
! ! ! ! Compute gridded area fractions:
! ! ! ! call par2grid

        ! Form divergence field * dt and store in ud temporarily:
        ! (normalize volg by cell volume)
        ud = volg / vcell - one

!         ! Write data
!         write(80,rec=istep) real(istep),real(ud)

        !-----------------------------------------
        ! Forward z cosine FFT:
        call dct(nx, nz, ud, ztrig, zfactors)
        ! Transpose array:
        do ix = 1, nx
            do kz = 0, nz
                phi(kz, ix) = ud(ix, kz)
            enddo
        enddo
        ! Forward x FFT:
        call forfft(nz+1, nx, phi, xtrig, xfactors)

        ! Invert Laplace's operator spectrally:
        phi = laplinv * phi
        ! phi = (spectral) velocity potential * dt

        !-----------------------------------------
        ! Compute x derivative spectrally:
        call deriv(nz+1, nx, hrkx, phi, uds)

        ! Reverse x FFT:
        call revfft(nz+1, nx, uds, xtrig, xfactors)
        ! Transpose array:
        do kz = 0, nz
            do ix = 1, nx
                ud(ix, kz) = uds(kz, ix)
            enddo
        enddo
        ! Reverse z cosine FFT:
        call dct(nx, nz, ud, ztrig, zfactors)

        !-----------------------------------------
        ! Compute z derivative spectrally:
        do kx = 1, nx
            do kz = 1, nz
                wds(kz, kx) = -rkz(kz) * phi(kz, kx)
            enddo
        enddo
        ! This makes wds a sine series in z

        ! Reverse x FFT:
        call revfft(nz, nx, wds, xtrig, xfactors)
        ! Transpose array:
        do kz = 1, nz
            do ix = 1, nx
                wka(ix, kz) = wds(kz, ix)
            enddo
        enddo
        ! Reverse z sine FFT:
        call dst(nx, nz, wka, ztrig, zfactors)

        ! Copy into wd with zero edge values:
        wd(:, 0) = zero
        do iz = 1, nz-1
            wd(:, iz) = wka(:, iz)
        enddo
        wd(:, nz) = zero

        !------------------------------------------------------------------
        ! Increment parcel positions usind (ud,wd) field:
        do n = 1, n_parcels

            call trilinear(parcels%position(n, i), ij, weight, ngp)

            do i = 1, ngp
                parcels%position(n, 1) = parcels%position(n, 1)             &
                                       + weight(i) * ud(ij(1, i), ij(2, i))

                parcels%position(n, 2) = parcels%position(n, 2)             &
                                       + weight(i) * wd(ij(1, i), ij(2, i))
            enddo

            call apply_periodic_bc(parcels%position(n, :))
        enddo
    end subroutine apply_divergent_flow

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! subroutine par2grid
!
! ! Reverse bi-linearly interpolates the area fraction (v)
! ! on all parcels to obtain the gridded area fraction (volg).
!
!  !Declarations:
! implicit none
!
!  !Local variables:
! double precision:: xx,px,pxc,zz,pz,pzc
! integer:: i,ix0,ix1,iz0,iz1
!
! !--------------------------------------------------
! volg = zero
!
! do i = 1,n
!   xx = dxi*(x(i)-xmin)
!   ix0 = min(1+int(xx),nx)
!   pxc = dble(ix0)-xx
!   px = one-pxc
!   ix1 = ixp(ix0)
!
!   zz = dzi*(z(i)-zmin)
!   iz0 = min(int(zz),nz-1)
!   pz = zz-dble(iz0)
!   pzc = one-pz
!   iz1 = iz0+1
!
!   volg(ix0,iz0) = volg(ix0,iz0)+v(i)*pxc*pzc
!   volg(ix0,iz1) = volg(ix0,iz1)+v(i)*pxc*pz
!   volg(ix1,iz0) = volg(ix1,iz0)+v(i)*px*pzc
!   volg(ix1,iz1) = volg(ix1,iz1)+v(i)*px*pz
! enddo
!
!  !Double edge values:
! volg(:, 0) = two*volg(:, 0)
! volg(:,nz) = two*volg(:,nz)
!
!  !Write out diagnostics:
! write(40,'(f7.2,2(1x,f15.12))') t,sqrt(sum((volg-one)**2)/dble(ncell)), &
!                                   maxval(abs(volg-one))
!
! return
! end subroutine par2grid

end module parcel_diverge
