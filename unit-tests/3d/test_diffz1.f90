! =============================================================================
!                       Test subroutine diffz1
!
!  This unit test checks the subroutine diffz1 using the
!  function:
!               cos(k * x) * sin(l * y) * cos(m * z)
!  where k = 2pi/L_x, l = 2pi/L_y and m = pi/L_z and where x, y and z all start
!  at 0 (one could start at -pi for x and y just as well).
!  The subroutine diffz1 should return
!               - m * cos(k * x) * sin(l * y) * sin(m * z)
! =============================================================================
program test_diffz1
    use unit_test
    use constants, only : zero, one, two, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils, only : init_fft
    use inversion_mod, only : diffz1
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), ds(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, m, prefactor

    ! doubling nx, ny, nz reduces the error by a factor 16
    nx = 64
    ny = 64
    nz = 64
    lower  = (/zero, zero, zero/)
    extent =  (/pi, pi, pi/)

    allocate(fs(0:nz, nx, ny))
    allocate(ds(0:nz, nx, ny))
    allocate(ref_sol(0:nz, nx, ny))

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

                fs(iz, ix, iy) = dcos(k * x) * dsin(l * y) * dcos(m * z)
                ref_sol(iz, ix, iy) = -m * dcos(k * x) * dsin(l * y) * dsin(m * z)
            enddo
        enddo
    enddo

    call init_fft

    call diffz1(fs, ds)

    error = maxval(dabs(ds - ref_sol))

    call print_result_dp('Test inversion (diffz1)', error, atol=1.0e-7)

    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

end program test_diffz1
