!===================================================================
!         Finds the parcels nearest every "small" parcel

!         Initialise data first using ranpar.f90
!===================================================================

program nearest

 !Import parameters and constants:
use constants

implicit none

 !Total number of parcels:
integer:: n

 !Parcel positions (x,z) and area fractions:
double precision:: x(nm),z(nm),v(nm)

 !Used for searching for possible parcel merger:
integer:: nppc(ncell),kc1(ncell),kc2(ncell) !ncell = nx*nz
integer:: loc(nm),node(nm),isma(nm/8),ibig(nm/8),nmerge
logical:: merge(nm)

 !Other variables:
double precision:: vmin,vmax,delx,delz,dsq,dsqmin,xm,zm
integer:: i,ic,i0,imin,k,m
integer:: ix,iz,ix0,iz0

!---------------------------------------------------------------------
! Read parcel positions and volume fractions (of grid cell volume):
open(80,file='ini_parcels.asc',status='old')
read(80,*) n
do i=1,n
  read(80,*) x(i),z(i),v(i)
enddo
close(80)

!---------------------------------------------------------------------
! Ask for merger criterion:
vmin=minval(v(1:n))
vmax=maxval(v(1:n))
write(*,*)
write(*,'(a,f9.7)') ' Minimum area fraction of any parcel = ',vmin
write(*,'(a,f9.7)') ' Maximum area fraction of any parcel = ',vmax

write(*,*)
write(*,*) ' Area fraction below which a parcel should be merged?'
read(*,*) vmin

!POMP$ INST BEGIN(foo1)
! These parcels are marked for merger:
merge(1:n)=(v(1:n) < vmin)
nmerge=0
nppc=0 !nppc(ic) will contain the number of parcels in grid cell ic

! Form list of small parcels:
do i=1,n
  if (merge(i)) then
    nmerge=nmerge+1
    isma(nmerge)=i
  else
    ix=int(dxi*(x(i)-xmin))
    iz=int(dzi*(z(i)-zmin))

    ! Cell index of parcel:
    ic=1+ix+nx*iz !This runs from 1 to ncell

    ! Accumulate number of parcels in this grid cell:
    nppc(ic)=nppc(ic)+1

    ! Store grid cell that this parcel is in:
    loc(i)=ic
  endif
enddo

! Find arrays kc1(ic) & kc2(ic) which indicate the parcels in grid cell ic
! through i = node(k), for k = kc1(ic),kc2(ic):
kc1(1)=1
do ic=1,ncell-1
  kc1(ic+1)=kc1(ic)+nppc(ic)
enddo

kc2=kc1-1
do i=1,n
  if (.not. merge(i)) then
    ic=loc(i)
    k=kc2(ic)+1
    node(k)=i
    kc2(ic)=k
  endif
enddo
!POMP$ INST END(foo1)

!POMP$ INST BEGIN(foo2)
!---------------------------------------------------------------------
! Now find the nearest grid point to each small parcel (to be merged)
! and search over the surrounding 8 grid cells for the closest parcel:
do m=1,nmerge
  i0=isma(m)
  xm=x(i0) !SB: store value for reuse, compiler may be picking up on this anyway
  zm=z(i0) !SB: store value for reuse, compiler may be picking up on this anyway
  ! Parcel i0 is small and should be merged; find closest other:
  ix0=mod(nint(dxi*(xm-xmin)),nx) ! ranges from 0 to nx-1
  iz0=nint(dzi*(zm-zmin))         ! ranges from 0 to nz
  ! Grid point (ix0,iz0) is closest to parcel i0

  ! Initialise minimum squared distance between parcels and parcel index:
  dsqmin=ellx*ellx
  imin=0

  ! Loop over 8 cells surrounding (ix0,iz0):
  do iz=max(0,iz0-1),min(nz-1,iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
    do ix=ix0-1,ix0
      ! Cell index (accounting for x periodicity):
      ic=1+mod(nx+ix,nx)+nx*iz
      ! Search parcels for closest:
      do k=kc1(ic),kc2(ic)
        i=node(k)
        delz=z(i)-zm
        dsq=delz*delz
        if(dsq<dsqmin) then !SB: perform initial test to reject parcel and avoid working out delx
           delx=x(i)-xm
           delx=delx-ellx*dble(int(delx*hlxi)) ! works across periodic edge
           dsq=dsq+delx*delx
           if (dsq < dsqmin) then
              dsqmin=dsq
              imin=i
           endif
        end if
      enddo
    enddo
  enddo
  ! Store the index of the parcel to be merged with:
  ibig(m)=imin !SB: only update here
enddo
!POMP$ INST END(foo2)

!---------------------------------------------------------------------
! Output data:
open(31,file='domain.asc',status='replace')
write(31,*) nx,nz
write(31,*) ellx,ellz
close(31)

open(31,file='merge.asc',status='replace')
do m=1,nmerge
  write(31,*) isma(m),ibig(m)
enddo
close(31)

write(*,*)
write(*,*) ' A list of the small & large parcels to be merged is in merge.asc'
write(*,*)

 !End main program
end program nearest
