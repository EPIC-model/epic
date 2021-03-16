!===================================================================
!   Randomly places parcels with variable volume in the domain
!===================================================================

program ranpar

 !Import parameters and constants:
use constants

implicit none

 !Total number of parcels:
integer:: n

 !Parcel positions (x,z) and area fractions:
double precision:: x(nm),z(nm),v(nm)

 !Other variables:
double precision:: vrange,fac
integer, dimension(:), allocatable :: seed
integer:: i,k,ngen

!---------------------------------------------------------------------
! Generate parcel positions and volumes:
write(*,*) ' Enter the average number of parcels per grid cell:'
read(*,*) fac
n=nint(fac*dble(ncell))

write(*,*) ' We assume any parcel area is equally probable.'
write(*,*) ' Enter A_max/A_min:'
read(*,*) vrange

write(*,*) ' Enter an integer seed for the random number generator:'
read(*,*) ngen
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=ngen
do i=1,ngen
  call random_seed(put=seed)
enddo

vrange=vrange-1.d0
do i=1,n
  call random_number(fac)
  x(i)=ellx*fac+xmin
  call random_number(fac)
  z(i)=ellz*fac+zmin
  call random_number(fac)
  v(i)=vrange*fac+1.d0
enddo

! Re-normalise so that sum(v) = ncell = nx*nz:
fac=dble(ncell)/sum(v(1:n))
v(1:n)=fac*v(1:n)

!---------------------------------------------------------------------
! Write parcel positions and area fractions (of grid cell volume):
open(80,file='ini_parcels.asc',status='replace')
write(80,*) n
do i=1,n
  write(80,*) x(i),z(i),v(i)
enddo
close(80)

write(*,*)
write(*,*) ' The initial parcels are in ini_parcels.asc'
write(*,*)

 !End main program
end program ranpar
