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
    use constants, only : pi, zero, one, two, f12, f14, f32
    use parcel_container
    use parcel_interpl, only : vol2grid
    use parameters, only : lower, extent, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : volg
    use timer
    implicit none

    integer :: i, j, k, l, kk, jj, ii
    double precision, parameter :: angle = f12 * pi
    double precision, parameter :: lam = 3.5d0 ! >= 3.5 --> 1 ellipse point outside domain
    double precision :: error

    nx = 4
    ny = 4
#ifdef ENABLE_3D
    nz = 4
    lower  = (/-f32, -f32, -f32/)
    extent = (/0.4d0, 0.4d0, 0.4d0/)
#else
    nz = 0
    lower  = (/-f32, -f32/)
    extent = (/0.4d0, 0.4d0/)
#endif

    call update_parameters

    allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

    call parcel_alloc(64)


    l = 1
    do k = 0, nz-1
        do j = 0, ny-1
            do i = 0, nx-1
                do kk = 1, 4, 2
                    do jj = 1, 4, 2
                        do ii = 1, 4, 2
                            parcels%position(l, 1) = lower(1) + i * dx(1) + f14 * dx(1) * ii
                            parcels%position(l, 2) = lower(2) + i * dx(2) + f14 * dx(2) * jj
#ifdef ENABLE_3D
                            parcels%position(l, 3) = lower(3) + k * dx(3) + f14 * dx(3) * kk
#endif
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    n_parcels = 64

    volg = zero

    parcels%volume = f14 * vcell

    ! b11
    parcels%B(:, 1) = lam * dcos(angle) ** 2 + one / lam * dsin(angle) ** 2

    ! b12
    parcels%B(:, 2) = f12 * (lam - one / lam) * dsin(two * angle)

    ! FIXME
    parcels%B(:, 3:5) = zero

    call vol2grid


    error = abs(sum(volg(0:nz, 0:ny-1, 0:nx-1)) - ngrid * vcell)

    call print_result_dp('Test free slip', error)

    call parcel_dealloc
    deallocate(volg)

end program test_free_slip
