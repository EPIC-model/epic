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
!  For this setting the theoretical vorticity tendency is
!  vtend(z, y, x, 1) = alpha*k*m**2 * (k^2+l^2)^(-1)*sin(k*x+l*y)*cos(k*x+l*y)
!  vtend(z, y, x, 2) = alpha*l*m**2 * (k^2+l^2)^(-1)*sin(k*x+l*y)*cos(k*x+l*y)
!  vtend(z, y, x, 3) = -alpha * m * sin(m*z) * cos(m*z)
! =============================================================================
program test_vtend
    use unit_test
    use constants, only : one, two, pi, f12, f34, three
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, vtend, field_default
    use inversion_utils, only : init_fft, fftxyp2s
    use inversion_mod, only : vor2vel, vor2vel_timer, vorticity_tendency, vtend_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: vtend_ref(:, :, :, :)
    integer                       :: ix, iy, iz, ik, il, im
    double precision              :: x, y, z, alpha, fk2l2, k, l, m
    double precision              :: cosmz, sinmz, sinkxly, coskxly

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('vtend', vtend_timer)

    nx = 128
    ny = 128
    nz = 128

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)


    allocate(vtend_ref(-1:nz+1, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    call init_fft

    do ik = 1, 3
        k = dble(ik)
        do il = 1, 3
            l = dble(il)
            do im = 1, 3, 2
                m = dble(im)

                alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
                fk2l2 = one / dble(k ** 2 + l ** 2)

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
                            velog(iz, iy, ix, 3) = cosmz * coskxly

                            ! vorticity
                            vortg(iz, iy, ix, 1) = alpha * velog(iz, iy, ix, 1)
                            vortg(iz, iy, ix, 2) = alpha * velog(iz, iy, ix, 2)
                            vortg(iz, iy, ix, 3) = alpha * velog(iz, iy, ix, 3)

                            ! reference solution
                            vtend_ref(iz, iy, ix, 1) = alpha * k * m ** 2 * fk2l2 * sinkxly * coskxly
                            vtend_ref(iz, iy, ix, 2) = alpha * l * m ** 2 * fk2l2 * sinkxly * coskxly
                            vtend_ref(iz, iy, ix, 3) = -alpha * m * sinmz * cosmz
                        enddo
                    enddo
                enddo

                call vorticity_tendency

                error = max(error, maxval(dabs(vtend_ref(0:nz, :, :, :) - vtend(0:nz, :, :, :))))
            enddo
        enddo
    enddo

    call print_result_dp('Test vorticity tendency', error, atol=5.0e-2)

    deallocate(vtend_ref)

end program test_vtend
