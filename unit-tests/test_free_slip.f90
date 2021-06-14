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
    use unit_test
    use constants, only : pi, zero, one, two
    use parcel_container
    use parcel_interpl, only : vol2grid
    use options, only : parcel
    use parameters, only : lower, extent, update_parameters, vcell, dx, nx, nz, ngrid
    use fields, only : volg
    implicit none

    integer :: i, j, k, jj, ii
    double precision, parameter :: angle = 0.5d0 * pi
    double precision, parameter :: lam = 3.5d0 ! >= 3.5 --> 1 ellipse point outside domain
    double precision :: error

    nx = 4
    nz = 4
    lower  = (/-1.5d0, -1.5d0/)
    extent = (/0.4d0, 0.4d0/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:nx-1))

    call parcel_alloc(64)


    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + 0.25d0 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + 0.25d0 * dx(2) * jj
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    n_parcels = 64

    volg = zero

    parcel%is_elliptic = .true.

    parcels%volume = 0.25d0 * vcell

    ! b11
    parcels%B(:, 1) = lam * dcos(angle) ** 2 + one / lam * dsin(angle) ** 2

    ! b12
    parcels%B(:, 2) = 0.5d0 * (lam - one / lam) * dsin(two * angle)


    call vol2grid


    error = abs(sum(volg(0:nz, 0:nx-1)) - ngrid * vcell)

    call print_result_dp('Test free slip', error)

    call parcel_dealloc
    deallocate(volg)

end program test_free_slip
