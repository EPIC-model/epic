! =============================================================================
!                           Test FFT module
!
!                      Tests the indexing in FFTs.
!                    (start index 0 vs start index 1)
! =============================================================================
program test_fft_1

    !Import FFT library:
    use stafft, only : initfft, dct, forfft
    use deriv1d

    implicit none

    !Parameters and constants:
    !---------------------------------------------------------------------
    !Generic double precision numerical constants:
    double precision,parameter:: zero=0.d0, one=1.d0
    double precision,parameter:: two=2.d0

    double precision,parameter:: f12=one/two

    !Number of x & z grid divisions:
    integer,parameter:: nx=40, nz=20

    !Maximum number of parcels:
    integer,parameter:: nm=64*nx*nz

    !Box width and height:
    double precision,parameter:: ellx=4.d0, ellz=2.d0

     !Domain half widths and edge values:
    double precision,parameter:: hlx=ellx/two

    !Grid lengths and their inverses:
    double precision,parameter:: dx=ellx/dble(nx)
    double precision,parameter:: dz=ellz/dble(nz)


    !---------------------------------------------------------------------


    double precision:: f0(nx,0:nz),f1(nx,0:nz)
    double precision:: f2(0:nx-1,0:nz)
    double precision:: s1(0:nz,nx),s2(0:nz,0:nx-1)

    double precision:: hrkx(nx),rkx(nx),rkz(nz)

    double precision:: xtrig(2*nx),ztrig(2*nz)
    integer:: xfactors(5),zfactors(5)

    double precision:: xg,zg,emax
    integer:: ix,iz,kx,kz

    !---------------------------------------------------------------------
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

    ! For z derivatives of a cosine in z function:
    rkz(nz)=zero

    !-----------------------------------------------------------------
    ! Define a field:
    do iz=0,nz
        zg=dz*dble(iz)
        do ix=1,nx
            xg=dx*dble(ix-1)-hlx
            f0(ix,iz)=(zg*(ellz-zg))**2/(one-f12*cos(xg+one))
        enddo
    enddo
    f1=f0
    f2=f0

    !-----------------------------------------
    ! Forward z cosine FFT:
    call dct(nx,nz,f1,ztrig,zfactors)
    ! Transpose array:
    do ix=1,nx
        do kz=0,nz
            s1(kz,ix)=f1(ix,kz)
        enddo
    enddo
    ! Forward x FFT:
    call forfft(nz+1,nx,s1,xtrig,xfactors)

    !-----------------------------------------
    ! Forward z cosine FFT:
    call dct(nx,nz,f2,ztrig,zfactors)
    ! Transpose array:
    do ix=0,nx-1
        do kz=0,nz
            s2(kz,ix)=f2(ix,kz)
        enddo
    enddo
    ! Forward x FFT:
    call forfft(nz+1,nx,s2,xtrig,xfactors)

    !-----------------------------------------
    ! Check difference:
    emax=zero
    do kx=1,nx
        do kz=0,nz
            emax=max(emax,abs(s1(kz,kx)-s2(kz,kx-1)))
        enddo
    enddo

    ! final check
    if (emax > 1.0e-15) then
        print '(a19, a17)', 'Test FFT indexing:', 'FAILED'
    else
        print '(a19, a17)', 'Test FFT indexing:', 'PASSED'
    endif

end program test_fft_1
