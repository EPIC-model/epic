! =============================================================================
!       Uses a divergent flow to push parcels toward or away from grid
!       points so as to bring the area fraction at each grid point back
!       to unity.
! =============================================================================
module parcel_diverge

    ! import FFT library:
    use stafft
    use deriv1d

    ! import parameters and constants:
    use constants

    implicit none

 !Total number of parcels:
integer:: n

 !Parcel positions (x,z) and area fractions:
double precision:: x(nm),z(nm),v(nm)

 !Rates of change (dx/dt,dz/dt):
double precision:: dxdt(nm),dzdt(nm)

 !Gridded area fraction:
double precision:: volg(nx,0:nz)

 !Wavenumbers and inverse Laplacian:
double precision:: hrkx(nx),rkx(nx),rkz(nz)
double precision:: laplinv(0:nz,nx)

 !Quantities needed in FFTs:
double precision:: xtrig(2*nx),ztrig(2*nz)
integer:: xfactors(5),zfactors(5)

 !Next horizontal grid point used in bilinear interpolation:
integer:: ixp(nx)

 !Time step (does not matter):
double precision,parameter:: dt=0.01d0

double precision:: t
integer:: istep,nsteps

!---------------------------------------------------------------------
! Initialise (read data, define fixed arrays):
call initialise

write(*,*) ' Enter the number of time steps used to relax the'
write(*,*) ' gridded area fraction back to unity:'
read(*,*) nsteps

! Evolve parcels:
do istep=1,nsteps
  call evolve(istep)
enddo

! Finalise:
call par2grid

! Write final data
write(80,rec=nsteps+1) real(nsteps+1),real(volg-one)

close(40)
close(80)

write(*,*)
write(*,*) ' The rms area fraction error is in volf.asc and the'
write(*,*) ' gridded area fraction error is in dv.r4; view with'
write(*,*)
write(*,'(a,i3,1x,i3,a)') '  dataview -ndim ',nx,nz+1,' dv.r4 -mod &'
write(*,*)

contains

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine initialise
! Initialises this program

implicit none

integer:: i,ix,kx,kz

!---------------------------------------------------------------------
! Read parcel positions and volume fractions (of grid cell volume):
open(80,file='ini_parcels.asc',status='old')
read(80,*) n
do i=1,n
  read(80,*) x(i),z(i),v(i)
enddo
close(80)

t=zero

!----------------------------------------------------------
! Next grid point indices used in bi-linear interpolation:
do ix=1,nx-1
  ixp(ix)=ix+1
enddo
ixp(nx)=1

!----------------------------------------------------------
! Set up FFTs:
call initfft(nx,xfactors,xtrig)
call initfft(nz,zfactors,ztrig)

! Define x wavenumbers:
call init_deriv(nx,ellx,hrkx)
rkx(1)=zero
do kx=1,nx/2-1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nx/2+1)=hrkx(nx)

! Define z wavenumbers:
call init_deriv(nz,ellz,rkz)

! Define spectral inverse Laplacian for inverting Poisson's equation:
do kx=1,nx
  do kz=1,nz
    laplinv(kz,kx)=-one/(rkx(kx)**2+rkz(kz)**2)
  enddo
enddo
! kz = 0:
do kx=2,nx
  laplinv(0,kx)=-one/rkx(kx)**2
enddo
! The zero wavenumber mode has no significance:
laplinv(0,1)=zero
! For z derivatives of a cosine in z function:
rkz(nz)=zero

!---------------------------------------------------------------
! File to contain rms area fraction - 1 & max abs value vs time:
open(40,file='volf.asc',status='replace')

! File to contain the gridded area fraction - 1 vs time:
open(80,file='dv.r4',form='unformatted',access='direct', &
                   status='replace',recl=4*(nx*(nz+1)+1))

return
end subroutine initialise

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine evolve(istep)
! Evolves the parcel positions (x,z) one time step (Euler method).

implicit none

double precision:: ud(nx,0:nz),wd(nx,0:nz),wka(nx,nz)
double precision:: phi(0:nz,nx),uds(0:nz,nx),wds(nz,nx)
double precision:: xx,px,pxc,zz,pz,pzc
integer:: istep,ix,iz,kx,kz
integer:: i,ix0,ix1,iz0,iz1

!-----------------------------------------------------------------
! Compute gridded area fractions:
call par2grid

! Form divergence field * dt and store in ud temporarily:
ud=volg-one

! Write data
write(80,rec=istep) real(istep),real(ud)

!-----------------------------------------
! Forward z cosine FFT:
call dct(nx,nz,ud,ztrig,zfactors)
! Transpose array:
do ix=1,nx
  do kz=0,nz
    phi(kz,ix)=ud(ix,kz)
  enddo
enddo
! Forward x FFT:
call forfft(nz+1,nx,phi,xtrig,xfactors)

! Invert Laplace's operator spectrally:
phi=laplinv*phi
! phi = (spectral) velocity potential * dt

!-----------------------------------------
! Compute x derivative spectrally:
call deriv(nz+1,nx,hrkx,phi,uds)

! Reverse x FFT:
call revfft(nz+1,nx,uds,xtrig,xfactors)
! Transpose array:
do kz=0,nz
  do ix=1,nx
    ud(ix,kz)=uds(kz,ix)
  enddo
enddo
! Reverse z cosine FFT:
call dct(nx,nz,ud,ztrig,zfactors)

!-----------------------------------------
! Compute z derivative spectrally:
do kx=1,nx
  do kz=1,nz
    wds(kz,kx)=-rkz(kz)*phi(kz,kx)
  enddo
enddo
! This makes wds a sine series in z

! Reverse x FFT:
call revfft(nz,nx,wds,xtrig,xfactors)
! Transpose array:
do kz=1,nz
  do ix=1,nx
    wka(ix,kz)=wds(kz,ix)
  enddo
enddo
! Reverse z sine FFT:
call dst(nx,nz,wka,ztrig,zfactors)

! Copy into wd with zero edge values:
wd(:,0)=zero
do iz=1,nz-1
  wd(:,iz)=wka(:,iz)
enddo
wd(:,nz)=zero

!------------------------------------------------------------------
! Increment parcel positions usind (ud,wd) field:
do i=1,n
  xx=dxi*(x(i)-xmin)
  ix0=min(1+int(xx),nx)
  pxc=dble(ix0)-xx
  px=one-pxc
  ix1=ixp(ix0)

  zz=dzi*(z(i)-zmin)
  iz0=min(int(zz),nz-1)
  pz=zz-dble(iz0)
  pzc=one-pz
  iz1=iz0+1

  x(i)=x(i)+pxc*(pzc*ud(ix0,iz0)+pz*ud(ix0,iz1))+ &
            +px*(pzc*ud(ix1,iz0)+pz*ud(ix1,iz1))

  z(i)=z(i)+pxc*(pzc*wd(ix0,iz0)+pz*wd(ix0,iz1))+ &
            +px*(pzc*wd(ix1,iz0)+pz*wd(ix1,iz1))
enddo

! Advance time:
t=t+dt

return
end subroutine evolve

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine par2grid

! Reverse bi-linearly interpolates the area fraction (v)
! on all parcels to obtain the gridded area fraction (volg).

 !Declarations:
implicit none

 !Local variables:
double precision:: xx,px,pxc,zz,pz,pzc
integer:: i,ix0,ix1,iz0,iz1

!--------------------------------------------------
volg=zero

do i=1,n
  xx=dxi*(x(i)-xmin)
  ix0=min(1+int(xx),nx)
  pxc=dble(ix0)-xx
  px=one-pxc
  ix1=ixp(ix0)

  zz=dzi*(z(i)-zmin)
  iz0=min(int(zz),nz-1)
  pz=zz-dble(iz0)
  pzc=one-pz
  iz1=iz0+1

  volg(ix0,iz0)=volg(ix0,iz0)+v(i)*pxc*pzc
  volg(ix0,iz1)=volg(ix0,iz1)+v(i)*pxc*pz
  volg(ix1,iz0)=volg(ix1,iz0)+v(i)*px*pzc
  volg(ix1,iz1)=volg(ix1,iz1)+v(i)*px*pz
enddo

 !Double edge values:
volg(:, 0)=two*volg(:, 0)
volg(:,nz)=two*volg(:,nz)

 !Write out diagnostics:
write(40,'(f7.2,2(1x,f15.12))') t,sqrt(sum((volg-one)**2)/dble(ncell)), &
                                  maxval(abs(volg-one))

return
end subroutine par2grid

end module parcel_diverge
