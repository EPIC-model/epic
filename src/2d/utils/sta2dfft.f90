module sta2dfft
    ! This module performs FFTs in two directions on two dimensional arrays using
    ! the stafft library module to actually compute the FFTs. If FFTs in one
    ! direction only are required use the stafft module directly. The module can
    ! compute any combination of sine, cosine and full FFTs in each direction.
    ! Along with the usual forwards (physical -> Fourier space) and reverse
    ! (Fourier space -> physical) routines there are also routines for computing
    ! the first derivatives in either direction.
    !
    ! The convention is that for each direction the array is dimensioned 1:nx or
    ! 1:ny for either the sine or full transforms. While the cosine transforms
    ! require the additional endpoint so 0:nx or 0:ny.
    !
    ! The routines contained in this module are:
    !
    ! init2dfft(nx, ny, lx, ly, xfactors, yfactors, xtrig, ytrig, kx, ky)
    !          This routine initialises all the arrays needed for further
    !          transforms. The integers nx and ny are the array dimensions. Then
    !          lx and ly are the domain lengths - these are needed for the correct
    !          scaling when computing derivatives. The arrays xfactors, yfactors,
    !          xtrig and ytrig are needed to perform the various FFTs by the stafft
    !          module (see there for further details. kx and ky are arrays to hold
    !          the wavenumbers associated with each mode in the domain, and are
    !          used in computing derivatives.
    !
    !          **If it is known at initialisation that no derivatives are required
    !            it is possible just to pass one for each of lx and ly, along with
    !            dummy arrays for kx and ky since these are only needed for
    !            computing the derviatives.**


    use constants, only : pi, zero
    use stafft

    implicit none

contains

    ! This subroutine performs the initialisation work for all subsequent
    ! transform and derivative routines.
    ! It calls the initfft() routine from the supproting 1d FFT module for
    ! transforms in both x and y directions.
    ! The routine then defines the two wavenumber arrays, one in each direction.
    subroutine init2dfft(nx, ny, lx, ly, xfactors, yfactors, xtrig, ytrig, kx, ky)
        integer,          intent(in)  :: nx, ny
        double precision, intent(in)  :: lx, ly
        integer,          intent(out) :: xfactors(5), yfactors(5)
        double precision, intent(out) :: xtrig(2 * nx), ytrig(2 * ny)
        double precision, intent(out) :: kx(nx), ky(ny)

        !Local declarations:
        double precision:: sc
        integer:: k

        !----------------------------------------------
        call initfft(nx, xfactors, xtrig)
        call initfft(ny, yfactors, ytrig)

        if (lx .ne. zero) then
            !Define x wavenumbers:
            sc = pi / lx
            do k = 1, nx
                kx(k) = sc * dble(k)
            enddo
        else
            !Catastrophic end to run if wave number definition fails:
            write(*, *) '**********************************************'
            write(*, *) ' Wavenumber array definition not possible.'
            write(*, *) ' Domain length in x equal to zero not allowed.'
            write(*, *) ' STOPPING...'
            write(*, *) '**********************************************'
            stop
        endif

        if (ly .ne. zero) then
            !Define y wavenumbers:
            sc = pi / ly
            do k = 1, ny
                ky(k) = sc * dble(k)
            enddo
        else
            !Catastrophic end to run if wave number definition fails:
            write(*, *) '**********************************************'
            write(*, *) ' Wavenumber array definition not possible.'
            write(*, *) ' Domain length in y equal to zero not allowed.'
            write(*, *) ' STOPPING...'
            write(*, *) '**********************************************'
            stop
        endif
    end subroutine

    !===================================================
    ! These routines are for arrays periodic in x and y.
    !===================================================

    ! Performs a physical -> spectral transform of a variable
    ! rvar(ny, nx) periodic in x and y, and returns the result
    ! (transposed) in svar(nx, ny).
    ! *** Note rvar is destroyed on return. ***
    subroutine ptospc(nx, ny, rvar, svar, xfactors, yfactors, xtrig, ytrig)
        integer,          intent(in)    :: nx, ny, xfactors(5), yfactors(5)
        double precision, intent(inout) :: xtrig(2 * nx), ytrig(2 * ny)
        double precision, intent(inout) :: rvar(ny, nx), svar(nx, ny)
        !Local declarations:
        integer:: kx, iy

        !Carry out a full x transform first:
        call forfft(ny, nx, rvar, xtrig, xfactors)

        !Transpose array:
        do kx = 1, nx
            do iy = 1, ny
                svar(kx, iy) = rvar(iy, kx)
            enddo
        enddo

        !Carry out a full y transform on transposed array:
        call forfft(nx, ny, svar, ytrig, yfactors)
    end subroutine


    !===================================================
    ! These routines are for arrays periodic in x only.
    !===================================================

    ! Performs a physical -> spectral transform of a variable
    ! rvar(0:ny, nx) periodic in x and represented by a cosine
    ! series in y and returns the result (transposed) in
    ! svar(nx, 0:ny).
    ! *** Note rvar is destroyed on return. ***
    subroutine ptospc_fc(nx, ny, rvar, svar, xfactors, yfactors, xtrig, ytrig)
        integer,          intent(in)    :: nx, ny, xfactors(5), yfactors(5)
        double precision, intent(inout) :: xtrig(2 * nx), ytrig(2 * ny)
        double precision, intent(inout) :: rvar(0:ny, nx)
        double precision, intent(out)   :: svar(nx, 0:ny)
        !Local declarations:
        integer:: kx, iy

        !Carry out a full x transform first:
        call forfft(ny + 1, nx, rvar, xtrig, xfactors)

        !Transpose array:
        do kx = 1, nx
            do iy = 0, ny
                svar(kx, iy) = rvar(iy, kx)
            enddo
        enddo

        !Carry out y cosine transform on transposed array:
        call dct(nx, ny, svar, ytrig, yfactors)
    end subroutine

end module
