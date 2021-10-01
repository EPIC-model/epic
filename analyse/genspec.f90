program genspec

use sta2dfft
use constants, only : pi, twopi

 !Declarations:
implicit none

 !Grid dimensions:
integer,parameter:: nx=64, ny=32

 !Width and height of the domain:
double precision,parameter:: ellx=1.d0, elly=1.d0

 !Array to contain data:
double precision:: pp(0:ny,0:nx-1)
 !Its Fourier transform:
double precision:: ss(0:nx-1,0:ny)
 !The spectrum:
double precision:: spec(0:max(nx,ny))

!x and y wavenumbers:
double precision:: rkx(0:nx-1),hrkx(nx),rky(ny)
 !Generic arrays needed for the FFTs:
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)
 !Wavenumber magnitude used for the spectrum:
integer:: kmag(0:nx-1,0:ny)

 !Other work variables:
double precision:: scx,rkxmax,scy,rkymax,delki
integer:: kxc,kmax,kx,ky,k

!---------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

 !Define x wavenumbers:
rkx(0)=0.d0
do kx=1,nx/2-1
  kxc=nx-kx
  rkx(kx )=hrkx(2*kx)
  rkx(kxc)=hrkx(2*kx)
enddo
rkx(nx/2)=hrkx(nx)

 !Initialise arrays for computing the spectrum:
scx=twopi/ellx
rkxmax=scx*dble(nx/2)
scy=pi/elly
rkymax=scy*dble(ny)
delki=1.d0/sqrt(scx**2+scy**2)
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do ky=1,ny
  do kx=0,nx-1
    kmag(kx,ky)=nint(sqrt(rkx(kx)**2+rky(ky)**2)*delki)
  enddo
enddo
do kx=0,nx-1
  kmag(kx,0)=nint(rkx(kx)*delki)
enddo

!---------------------------------------------------------------------
 !Read data into array pp:



!---------------------------------------------------------------------
 !Compute spectrum:

 !Transform data in pp to spectral space:
call ptospc_fc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)

do k=0,kmax
  spec(k)=0.d0
enddo

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+0.25d0*ss(0,0)**2

 !y-independent mode:
do kx=1,nx-1
  k=kmag(kx,0)
  spec(k)=spec(k)+0.5d0*ss(kx,0)**2
enddo

 !x-independent mode:
do ky=1,ny
  k=kmag(0,ky)
  spec(k)=spec(k)+0.5d0*ss(0,ky)**2
enddo

 !All other modes:
do ky=1,ny
  do kx=1,nx-1
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

!---------------------------------------------------------------------
 !Write spectrum contained in spec(k):




end program genspec
