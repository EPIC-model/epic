! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz using the
!  function:
!               cos(k * x) * sin(l * y) * sin(m * z)
!  where k = 2pi/L_x, l = 2pi/L_y and m = pi/L_z and where x, y and z all start
!  at 0 (one could start at -pi for x and y just as well).
!  The subroutine diffz should return
!               m * cos(k * x) * sin(l * y) * cos(m * z)
! =============================================================================
program test_diffz
    use unit_test
    use constants, only : zero, one, two, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils, only : init_inversion, field_combine_physical, field_decompose_physical
    use inversion_mod, only : diffz
    implicit none

    double precision              :: error
    double precision, allocatable :: fp(:, :, :), dp(:, :, :), fs(:, :, :), ds(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, m, prefactor

    nx = 32
    ny = 64
    nz = 128
    lower  = (/zero, zero, zero/)
    extent =  (/pi, twopi, two * twopi/)

    allocate(fp(0:nz, ny, nx))
    allocate(dp(0:nz, ny, nx))
    allocate(fs(0:nz, nx, ny))
    allocate(ds(0:nz, nx, ny))
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
                ref_sol(iz, iy, ix) = m * dcos(k * x) * dsin(l * y) * dcos(m * z)

            enddo
        enddo
    enddo

    call init_inversion

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)

    error = maxval(dabs(dp - ref_sol))

    call print_result_dp('Test inversion (diffz)', error, atol=1.5e-13)

    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

end program test_diffz
