!===================================================================
!  Uses a divergent flow to push parcels toward or away from grid
!  points so as to bring the area fraction at each grid point back
!  to unity.

!  Here tri-diagonal solves are used in the vertical direction.

!           Initialise data first using inipar.f90
!===================================================================

program diverge

 !Import FFT library:
use stafft
use deriv1d

 !Import parameters and constants:
use constants

implicit none

 !Total number of parcels:
integer:: n

 !Parcel positions (x,z) and area fractions:
double precision:: x(nm),z(nm),v(nm)

 !Rates of change (dx/dt,dz/dt):
double precision:: dxdt(nm),dzdt(nm)

 !Gridded area fraction:
double precision:: volg(0:nz,nx)

!Tri-diagonal arrays:
double precision:: ap(nx),apb(nx)
double precision:: etdv(0:nz,nx),htdv(0:nz,nx)
double precision:: etd1(nz-1),htd1(nz-1)
double precision:: etda(nz),htda(nz)

 !Wavenumbers::
double precision:: hrkx(nx),rkx(nx)

 !Quantities needed in FFTs:
double precision:: xtrig(2*nx)
integer:: xfactors(5)

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

double precision:: a0(nx),a0b(nx),ksq(nx)
double precision,parameter:: dzisq=dzi**2
integer,parameter:: nwx=nx/2
integer:: i,ix,kx,iz

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

! Define x wavenumbers:
call init_deriv(nx,ellx,hrkx)
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)

 !Squared wavenumber array (used in tridiagonal solve):
ksq=rkx**2

!-----------------------------------------------------------------------
! Fixed coefficients used in the tridiagonal problems:
a0=-two*dzisq-f56*ksq
a0b=-dzisq-f13*ksq
ap=dzisq-f112*ksq
apb=dzisq-f16*ksq

!-----------------------------------------------------------------------
! Tridiagonal arrays for inversion of Poisson's equation:
htdv(:,1)=zero
etdv(:,1)=zero
do kx=2,nx
  htdv(0,kx)=one/a0b(kx)
  etdv(0,kx)=-apb(kx)*htdv(0,kx)
  do iz=1,nz-1
    htdv(iz,kx)=one/(a0(kx)+ap(kx)*etdv(iz-1,kx))
    etdv(iz,kx)=-ap(kx)*htdv(iz,kx)
  enddo
  htdv(nz,kx)=one/(a0b(kx)+apb(kx)*etdv(nz-1,kx))
enddo

! Tridiagonal arrays for the compact difference calculation of d/dz
! for fields f for which df/dz = 0 at the boundaries:
htd1(1)=one/f23
etd1(1)=-f16*htd1(1)
do iz=2,nz-2
  htd1(iz)=one/(f23+f16*etd1(iz-1))
  etd1(iz)=-f16*htd1(iz)
enddo
htd1(nz-1)=one/(f23+f16*etd1(nz-2))

! Tridiagonal arrays used for integrating in z (see vertint):
htda(1)=one/f76
etda(1)=-f16*htda(1)
do iz=2,nz-1
  htda(iz)=one/(one+f16*etda(iz-1))
  etda(iz)=-f16*htda(iz)
enddo
htda(nz)=one/(f76+f16*etda(nz-2))

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

double precision:: phi(0:nz,nx),ud(0:nz,nx),wd(0:nz,nx)
double precision:: wbar(0:nz)
double precision:: xx,px,pxc,zz,pz,pzc
integer:: istep,ix,iz,kx
integer:: i,ix0,ix1,iz0,iz1
!DGD
double precision,parameter:: pi=acos(-1.d0)
double precision,parameter:: k=two*pi/ellx,m=pi/ellz
double precision:: phie(0:nz,nx),ude(0:nz,nx),wde(0:nz,nx)
!DGD

!-----------------------------------------------------------------
! Compute gridded area fractions:
call par2grid

! Form divergence field * dt and store in phi temporarily:
phi=volg-one

! Write data
write(80,rec=istep) real(istep),real(phi)

!DGD
!px=one
px=one
do ix=1,nx
  xx=k*(dx*dble(ix-1)-hlx)
  do iz=0,nz
    zz=m*(zmin+dz*dble(iz))
    phie(iz,ix)=sin(xx)*sin(zz)
    ude(iz,ix)=k*cos(xx)*sin(zz)
    wde(iz,ix)=m*(px+sin(xx))*cos(zz)
    phi(iz,ix)=-(m**2*px+(k**2+m**2)*sin(xx))*sin(zz)
  enddo
enddo
!DGD

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
wd(:,1)=wd(:,1)+wbar

! Reverse x FFT:
call revfft(nz+1,nx,wd,xtrig,xfactors)

!DGD
open(77,file='nduw.r4',form='unformatted',access='direct', &
                     status='replace',recl=4*(nx*(nz+1)+1))
write(77,rec=1) real(0.),real(ud-ude)
write(77,rec=2) real(0.),real(wd-wde)
close(77)

call revfft(nz+1,nx,phi,xtrig,xfactors)
open(77,file='ndphi.r4',form='unformatted',access='direct', &
                      status='replace',recl=4*(nx*(nz+1)+1))
write(77,rec=1) real(0.),real(phi-phie)
close(77)

open(77,file='phi.r4',form='unformatted',access='direct', &
                    status='replace',recl=4*(nx*(nz+1)+1))
write(77,rec=1) real(0.),real(phi)
close(77)

ud=(ud-ude)**2
xx=sqrt((f12*sum(ud(0,:)+ud(nz,:))+sum(ud(1:nz-1,:)))/dble(ncell))
wd=(wd-wde)**2
zz=sqrt(sum(wd(1:nz-1,:))/dble(ncell))
write(*,'(a,i4,2(a,f18.16))') &
     ' nx = ',nx,'   rms ud error = ',xx,'   rms wd error = ',zz
stop
!DGD

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

  x(i)=x(i)+pxc*(pzc*ud(iz0,ix0)+pz*ud(iz1,ix0))+ &
            +px*(pzc*ud(iz0,ix1)+pz*ud(iz1,ix1))

  z(i)=z(i)+pxc*(pzc*wd(iz0,ix0)+pz*wd(iz1,ix0))+ &
            +px*(pzc*wd(iz0,ix1)+pz*wd(iz1,ix1))
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

  volg(iz0,ix0)=volg(iz0,ix0)+v(i)*pxc*pzc
  volg(iz1,ix0)=volg(iz1,ix0)+v(i)*pxc*pz
  volg(iz0,ix1)=volg(iz0,ix1)+v(i)*px*pzc
  volg(iz1,ix1)=volg(iz1,ix1)+v(i)*px*pz
enddo

 !Double edge values:
do ix0=1,nx
  volg(0 ,ix0)=two*volg(0 ,ix0)
  volg(nz,ix0)=two*volg(nz,ix0)
enddo

 !Write out diagnostics:
write(40,'(f7.2,2(1x,f15.12))') t,sqrt(sum((volg-one)**2)/dble(ncell)), &
                                  maxval(abs(volg-one))

return
end subroutine par2grid

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine lapinv1(fs)
 !Inverts Laplace's operator on fs in semi-spectral space.
 !Here dfs/dz = 0 on the z boundaries.
 !Uses 4th-order compact differencing
 !*** Overwrites fs ***

 !Declarations:
implicit none

 !Passed variables:
double precision:: fs(0:nz,nx)

 !Local variables:
double precision:: rs(0:nz,nx)
integer:: kx,iz

fs(:,1)=zero
do kx=2,nx
  rs(0,kx)=f13*fs(0,kx)+f16*fs(1,kx)
  do iz=1,nz-1
    rs(iz,kx)=f112*(fs(iz-1,kx)+fs(iz+1,kx))+f56*fs(iz,kx)
  enddo
  rs(nz,kx)=f13*fs(nz,kx)+f16*fs(nz-1,kx)

  fs(0,kx)=rs(0,kx)*htdv(0,kx)
  do iz=1,nz-1
    fs(iz,kx)=(rs(iz,kx)-ap(kx)*fs(iz-1,kx))*htdv(iz,kx)
  enddo
  fs(nz,kx)=(rs(nz,kx)-apb(kx)*fs(nz-1,kx))*htdv(nz,kx)

  do iz=nz-1,0,-1
    fs(iz,kx)=etdv(iz,kx)*fs(iz+1,kx)+fs(iz,kx)
  enddo
enddo

return
end subroutine

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine diffz1(fs,ds)
 !Calculates df/dz for a field f which has df/dz = 0 at the boundaries
 !using 4th-order compact differencing.  Here fs = f and ds = df/dz.

 !Declarations:
implicit none

 !Passed variables:
double precision:: fs(0:nz,nx),ds(0:nz,nx)

 !Local variables:
integer:: ix,iz

do ix=1,nx
  do iz=1,nz-1
    ds(iz,ix)=(fs(iz+1,ix)-fs(iz-1,ix))*hdzi
  enddo

  ds(0,ix)=zero
  ds(1,ix)=ds(1,ix)*htd1(1)
  do iz=2,nz-1
    ds(iz,ix)=(ds(iz,ix)-f16*ds(iz-1,ix))*htd1(iz)
  enddo
  ds(nz,ix)=zero

  do iz=nz-2,1,-1
    ds(iz,ix)=etd1(iz)*ds(iz+1,ix)+ds(iz,ix)
  enddo
enddo

return
end subroutine

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine vertint(ds,fs)
 !Finds f by integrating df/dz = d, ensuring f = 0 at the boundaries
 !using 4th-order compact differencing.  Here ds = df/dz and fs = f.

 !Declarations:
implicit none

 !Passed variables:
double precision:: ds(0:nz),fs(0:nz)

 !Local variables:
double precision:: es(nz),esum
integer:: iz

!-------------------------------------------
 !First interpolate ds to a half grid as es:
do iz=1,nz
  es(iz)=f23*(ds(iz-1)+ds(iz))
enddo

es(1)=es(1)*htda(1)
do iz=2,nz
  es(iz)=(es(iz)-f16*es(iz-1))*htda(iz)
enddo

do iz=nz-1,1,-1
  es(iz)=etda(iz)*es(iz+1)+es(iz)
enddo

!-------------------------------------------
 !Next adjust es to ensure f(nz) = 0:
esum=(f1112*(es(1)+es(nz))+sum(es(2:nz-1)))/(dble(nz)-f16)
es=es-esum

 !Integrate:
fs(0)=zero
fs(1)=dz24*(23.d0*es(1)+es(2))
do iz=2,nz-1
  fs(iz)=fs(iz-1)+dz24*(es(iz-1)+22.d0*es(iz)+es(iz+1))
enddo
fs(nz)=zero

return
end subroutine

 !End main program
end program diverge
