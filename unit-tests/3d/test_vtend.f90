! =============================================================================
!                       Test the vorticity tendency
!
!  This unit test checks the calculation of the vorticity tendency using the
!  Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
!  The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
!  We use k = l = 2 and m = 1. For this setting the theoertical vorticity
!  tendency is
!        omega*grad(u) = (3/8)sin(2 theta) (ii + jj) - (3/2)sin(2 z) kk
!  where ii, jj and kk are the unit vectors and theta = 2(x+y)
! =============================================================================
program test_vtend
    use unit_test
    use constants, only : one, two, pi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, vtend, tbuoyg, field_default
    use inversion_utils, only : init_fft, fftxyp2s
    use inversion_mod, only : vor2vel, vor2vel_timer, vorticity_tendency, vtend_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: vtend_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, m, alpha, fk2l2
    double precision              :: cosmz, sinmz, sinkxly, coskxly

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('vtend', vtend_timer)

    nx = 64
    ny = 64
    nz = 64

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    k = two
    l = two
    m = one
    alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
    fk2l2 = alpha / dble(k ** 2 + l ** 2)

    allocate(vtend_ref(-1:nz+1, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = -1, nz+1
                z = lower(3) + iz * dx(3)

                cosmz = dcos(m * z)
                sinmz = dsin(m * z)
                sinkxly = dsin(k * x + l * y)
                coskxly = dcos(k * x + l * y)

                ! velocity
                velog(iz, iy, ix, 1) = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                velog(iz, iy, ix, 2) = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                velog(iz, iy, ix, 3) = alpha * cosmz * coskxly

                ! vorticity
                vortg(iz, iy, ix, 1) = alpha * velog(iz, iy, ix, 1)
                vortg(iz, iy, ix, 1) = alpha * velog(iz, iy, ix, 2)
                vortg(iz, iy, ix, 1) = alpha * velog(iz, iy, ix, 3)

                ! reference solution
                vtend_ref(iz, iy, ix, 1) = 3.0d0 / 8.0d0 * sinkxly
                vtend_ref(iz, iy, ix, 2) = 3.0d0 / 8.0d0 * sinkxly
                vtend_ref(iz, iy, ix, 3) = - 3.0d0 / 2.0d0 * dsin(2.0d0 * z)
            enddo
        enddo
    enddo

    call init_fft

    call vorticity_tendency(vortg, velog, tbuoyg, vtend)

    error = maxval(dabs(vtend_ref(0:nz, :, :, :) - vtend(0:nz, :, :, :)))

    print *, error

    call print_result_dp('Test vorticity tendency', error, atol=4.0e-7)

    deallocate(vtend_ref)

end program test_vtend
