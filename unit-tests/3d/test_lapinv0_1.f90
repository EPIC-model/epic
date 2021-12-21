! =============================================================================
!                       Test subroutine lapinv0
!
!  This unit test checks the subroutine lapinv0 using the the
!  function:
!              cos(k * x) * sin(l * y) * sin(m * z)
!  where k = 2pi/L_x, l = 2pi/L_y and m = pi/L_z and where x, y and z all start
!  at 0 (one could start at -pi for x and y just as well).
!  The subroutine lapinv0 should return the same function multiplied by
!  -1/(k^2 + l^2 + m^2).
! =============================================================================
program test_lapinv0_1
    use unit_test
    use constants, only : zero, one, two, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils, only : init_fft, fftxyp2s, fftxys2p
    use inversion_mod, only : lapinv0
    implicit none

    double precision              :: error
    double precision, allocatable :: fp(:, :, :), fs(:, :, :), ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, m, prefactor

    nx = 32
    ny = 64
    nz = 128
    lower  = (/zero, zero, zero/)
    extent =  (/pi, twopi, two * twopi/)

    allocate(fs(0:nz, nx, ny))
    allocate(fp(0:nz, ny, nx))
    allocate(ref_sol(0:nz, ny, nx))

    call update_parameters

    k = twopi / extent(1)
    l = twopi / extent(2)
    m =    pi / extent(3)

    prefactor = - one / (k ** 2 + l ** 2 + m ** 2)

    do ix = 1, nx
        x = lower(1) + (ix - 1) * dx(1)
        do iy = 1, ny
            y = lower(2) + (iy - 1) * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                fp(iz, iy, ix) = dcos(k * x) * dsin(l * y) * dsin(m * z)

                ref_sol(iz, iy, ix) = prefactor * fp(iz, iy, ix)
            enddo
        enddo
    enddo

    call init_fft

    call fftxyp2s(fp, fs)

    call lapinv0(fs)

    call fftxys2p(fs, fp)

    error = maxval(dabs(fp - ref_sol))

    call print_result_dp('Test inversion (lapinv0)', error, atol=4.0e-12)

    deallocate(fs)
    deallocate(fp)
    deallocate(ref_sol)

end program test_lapinv0_1
