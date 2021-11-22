! =============================================================================
!                       Test the inversion
!
!         This unit test checks the inversion of the vorticity.
! =============================================================================
program test_vor2vel
    use unit_test
    use constants, only : pi, zero, one, two, twopi, six, three, four, f12, f13, f14, f15
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, velgradg, field_alloc
    use inversion_utils, only : init_fft, fftxys2p
    use inversion_mod, only : vor2vel
    use timer
    implicit none

    double precision              :: error = zero
    double precision, parameter   :: k = one                 &
                                   , l = two                 &
                                   , kmag2 = k ** 2 + l ** 2 &
                                   , alpha = -one            &
                                   , beta = two              &
                                   , lz = twopi              &
                                   , p1 = 0.2d0 * lz         &
                                   , p2 = 0.7d0 * lz         &
                                   , p3 = lz                 &
                                   , q1 = 0.3d0 * lz         &
                                   , q2 = 0.8d0 * lz         &
                                   , q3 = lz
    double precision, allocatable :: svelog_ref(:, :, :, :), svortg(:, :, :, :), &
                                     vor_copy(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: Ahat, Bhat, xi, eta, zeta, pp, qq, ps, qs
    double precision              :: z, z2, z3, gam, C, zq1, zq2, zq3, zp1, zp2, zp3
    double precision              :: us, vs, ws, ppp, qqq, x, y

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-pi, -pi, zero/)
    extent =  (/lz, lz, lz/)

!     allocate(svelog(-1:nz+1, 0:nx-1, 0:ny-1, 3))
    allocate(svelog_ref(-1:nz+1, 0:nx-1, 0:ny-1, 3))
    allocate(svortg(0:nz, 0:nx-1, 0:ny-1, 3))
    allocate(vor_copy(0:nz, 0:nx-1, 0:ny-1, 3))

    call update_parameters

    call field_alloc

    pp = p1 * p2 + p1 * p3 + p2 * p3
    qq = q1 * q2 + q1 * q3 + q2 * q3
    ps = p1 + p2 + p3
    qs = q1 + q2 + q3

    ppp = p1 * p2 * p3
    qqq = q1 * q2 * q3

    gam = 0

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            gam = zero
            do iz = 0, nz
                z = iz * dx(3)
                z2 = z ** 2
                z3 = z2 * z

                gam = gam  - k * alpha * z2 * (f15 * z3 - f14 * ps * z2 + f13 * pp * z - f12 * ppp) &
                           - l * beta  * z2 * (f15 * z3 - f14 * qs * z2 + f13 * qq * z - f12 * qqq)
            enddo
            gam = gam / (nz+1)

            do iz = 0, nz
                z = iz * dx(3)
                z2 = z ** 2
                z3 = z2 * z

                zq1 = z - q1
                zq2 = z - q2
                zq3 = z - q3

                zp1 = z - p1
                zp2 = z - p2
                zp3 = z - p3

                Ahat = alpha * z * zp1 * zp2 * zp3 * dcos(k * x) * dcos(l * y)
                Bhat = beta  * z * zq1 * zq2 * zq3 * dcos(k * x) * dcos(l * y)

                C = gam - k * alpha * z2 * (f15 * z3 - f14 * ps * z2 + f13 * pp * z - f12 * ppp) &
                        - l * beta  * z2 * (f15 * z3 - f14 * qs * z2 + f13 * qq * z - f12 * qqq)

                xi  = two * alpha * (six * z2 - three * ps * z + pp) - kmag2 * Ahat
                eta = two * beta  * (six * z2 - three * qs * z + qq) - kmag2 * Bhat

                zeta = - k * alpha * (four * z3 - three * ps * z2 + two * pp * z - ppp) &
                       - l * beta  * (four * z3 - three * qs * z2 + two * qq * z - qqq) &
                       - kmag2 * C

                svortg(iz, ix, iy, 1) = xi
                svortg(iz, ix, iy, 2) = eta
                svortg(iz, ix, iy, 3) = zeta

                us =    beta * (z * zq1 * zq2 + z * zq3 * zq2 + zq1 * zq3 * zq2 + z * zq1 * zq3) + l * C
                vs = - alpha * (z * zp1 * zp2 + z * zp3 * zp2 + zp1 * zp3 * zp2 + z * zp1 * zp3) - k * C
                ws = l * Ahat - k * Bhat ! imaginary

                svelog_ref(iz, ix, iy, 1) = us
                svelog_ref(iz, ix, iy, 2) = vs
                svelog_ref(iz, ix, iy, 3) = ws
            enddo
        enddo
    enddo

    call init_fft

    vor_copy = svortg

    ! convert vortg to physical space
    call fftxys2p(svortg(0:nz, :, :, 1), vortg(0:nz, :, :, 1))
    call fftxys2p(svortg(0:nz, :, :, 2), vortg(0:nz, :, :, 2))
    call fftxys2p(svortg(0:nz, :, :, 3), vortg(0:nz, :, :, 3))

    call vor2vel(vortg, velog, velgradg)


!     error = maxval(dabs(svelog(0:nz, :, :, :) - svelog_ref(0:nz, :, :, :)))

!     do iz = 0, nz
!         z = iz * dx(3)
!         print *, z, vor_copy(iz, 15, 15, :), svelog_ref(iz, 15, 15, :), svelog(iz, 15, 15, :)
!     enddo
!     print *, error

!     call print_result_dp('Test inversion (vor2vel)', error)

!     deallocate(svelog)

end program test_vor2vel
