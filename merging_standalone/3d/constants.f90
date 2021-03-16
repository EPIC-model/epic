module constants
!==================================================================
!    This module contains all the non-modifiable parameters and
!    all quantities which never change throughout a simulation.
!==================================================================

 !Import parameters:
use parameters

integer,parameter:: nxm1=nx-1, nym1=ny-1, nzm1=nz-1

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0, one=1.d0
double precision,parameter:: two=2.d0, three=3.d0
double precision,parameter:: four=4.d0, six=6.d0

double precision,parameter:: f12=one/two, f13=one/three, f14=one/four
double precision,parameter:: f16=one/six, f56=5.d0/six, f112=one/12.d0

double precision,parameter:: small=1.d-12, oms=one-small

!---------------------------------------------------------------------
 !Domain half widths and edge values:
double precision,parameter:: hlx=ellx/two, hlxi=one/hlx, xmin=-hlx
double precision,parameter:: hly=elly/two, hlyi=one/hly, ymin=-hly
double precision,parameter:: zmin=-ellz/two, zmax=zmin+ellz

 !Grid lengths and their inverses:
double precision,parameter:: dx=ellx/dble(nx), dxi=dble(nx)/(ellx*oms)
double precision,parameter:: dy=elly/dble(ny), dyi=dble(ny)/(elly*oms)
double precision,parameter:: dz=ellz/dble(nz), dzi=dble(nz)/(ellz*oms)

 !Grid cell volume:
double precision,parameter:: vcell=dx*dy*dz

!---------------------------------------------------------------------
 !Number of grid cells:
integer,parameter:: ncell=nx*ny*(nz+1)

end module
