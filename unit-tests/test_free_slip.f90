! =============================================================================
!                       Test free slip boundary condition
!
!       This unit test checks the free slip boundary condition which is
!       implied in vertical direction. Ellipses with aspect ratio 3.5
!       are place such that the parcels at the bottom and ad the top have
!       one ellipse point reaching out the physical domain. We use 4
!       parcels per cell. After interpolating the parcel volume
!       (volume = 0.25 * vcell per parcel), the gridded volume is
!       ngrid * vcell where ngrid is the number of grid points.
! =============================================================================
program test_free_slip
    use constants, only : pi
    use parcel_container
    use interpolation, only : par2grid
    use options, only : parcel_info, grid
    use parameters, only : extent, lower, update_parameters, vcell, dx, nx, nz, ngrid
    implicit none

    double precision :: volg(-1:5, 0:3, 1)
    integer :: i, j, k, jj, ii
    double precision, parameter :: angle = 0.5 * pi
    double precision, parameter :: lam = 3.5 ! >= 3.5 --> 1 ellipse point outside domain
    double precision :: error

    grid = (/5, 5/)

    extent =  (/0.4, 0.4/)
    lower = (/-0.2, -0.2/)

    call update_parameters()

    call parcel_alloc(64)


    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + 0.25 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + 0.25 * dx(2) * jj
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    n_parcels = 64

    volg = 0.0

    parcel_info%is_elliptic = .true.

    parcels%volume = 0.25 * vcell

    ! b11
    parcels%B(:, 1) = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2

    ! b12
    parcels%B(:, 2) = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)


    call par2grid(parcels, parcels%volume, volg)


    error = abs(sum(volg(0:nz, 0:nx-1, :)) - ngrid * vcell)

    if (error > 1.0e-15) then
        print '(a16, a20)', 'Test free slip:', 'FAILED'
    else
        print '(a16, a20)', 'Test free slip:', 'PASSED'
    endif

end program test_free_slip
