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
program test_free_slip_3d
    use unit_test
    use constants, only : pi                      &
                        , zero , one, two, three  &
                        , f12, f14, f18, f23
    use parcel_container
    use parcel_ellipsoid, only : get_abc, get_B33
    use parcel_interpl, only : vol2grid
    use parameters, only : lower, extent, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : volg
    use timer
    implicit none

    integer :: i, j, k, l, kk, jj, ii
    double precision :: error, a2, b2, c2, abc, abc23

    nx = 2
    ny = 2
    nz = 2
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

    call parcel_alloc(2 * 8 * nx * ny * nz)

    l = 1
    do k = 0, nz-1
        do j = 0, ny-1
            do i = 0, nx-1
                do kk = 1, 4, 2
                    do jj = 1, 4, 2
                        do ii = 1, 4, 2
                            parcels%position(l, 1) = lower(1) + i * dx(1) + f14 * dx(1) * ii
                            parcels%position(l, 2) = lower(2) + j * dx(2) + f14 * dx(2) * jj
                            parcels%position(l, 3) = lower(3) + k * dx(3) + f14 * dx(3) * kk
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    n_parcels = l - 1

    volg = zero

    parcels%volume = f18 * vcell

    abc = get_abc(parcels%volume(1))
    abc23 = abc ** f23
    a2 = 1.5d0 * abc23
    b2 = abc23 / dsqrt(1.5d0)
    c2 = abc23 / dsqrt(1.5d0)

    parcels%B(:, 1) = c2
    parcels%B(:, 2) = zero
    parcels%B(:, 3) = zero
    parcels%B(:, 4) = b2
    parcels%B(:, 5) = zero
!     print *, get_B33(parcels%B(1, :), parcels%volume(1))
!     stop


!     do l = 1, n_parcels
!         print *, parcels%position(l, :), parcels%B(l, :), get_B33(parcels%B(l, :), parcels%volume(l))
!     enddo

!     stop

    call vol2grid

    print *, nx * ny * nz, ngrid, vcell * nx * ny * nz, sum(volg(0:nz, 0:ny-1, 0:nx-1))

    error = abs(sum(volg(0:nz, 0:ny-1, 0:nx-1)) - ngrid * vcell)

    print *, error

    call print_result_dp('Test free slip', error)

    call parcel_dealloc
    deallocate(volg)

end program test_free_slip_3d
