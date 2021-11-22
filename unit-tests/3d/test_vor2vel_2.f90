! =============================================================================
!                       Test the vorticity inversion
!
!  This unit test checks the calculation of the velocity using the
!  ABC flow:
!               u = A * sin(z) + C * cos(y)
!               v = B * sin(x) + A * cos(z)
!               w = C * sin(y) + B * cos(x)
!  The vorticity of this flow is
!               xi = A * sin(z) + C * cos(y)
!              eta = B * sin(x) + A * cos(z)
!             zeta = C * sin(y) + B * cos(x)
!  hence, the cross-product of vorticity and velocity (i.e. Lamb vector)
!  is zero.
!
!  Reference:
!  16 November 2021
!  https://en.wikipedia.org/wiki/Arnold%E2%80%93Beltrami%E2%80%93Childress_flow
! =============================================================================
program test_vor2vel_2
    use unit_test
    use constants, only : zero, two, four, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, velgradg, field_alloc
    use inversion_utils, only : init_fft
    use inversion_mod, only : vor2vel, vor2vel_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: velog_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, A, B, C

    call register_timer('vorticity', vor2vel_timer)


    nx = 32
    ny = 64
    nz = 128
    lower  = (/-pi, -twopi, -two * twopi/)
    extent =  (/twopi, two * twopi, four * twopi/)

    A = 0.5d0
    B = 1.2d0
    C = 0.9d0

    allocate(velog_ref(0:nz, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_alloc

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                ! vorticity
                vortg(iz, iy, ix, 1) = A * dsin(z) + C * dcos(y)
                vortg(iz, iy, ix, 2) = B * dsin(x) + A * dcos(z)
                vortg(iz, iy, ix, 3) = C * dsin(y) + B * dcos(x)

                ! velocity (reference solution)
                velog_ref(iz, iy, ix, 1) = vortg(iz, iy, ix, 1)
                velog_ref(iz, iy, ix, 2) = vortg(iz, iy, ix, 2)
                velog_ref(iz, iy, ix, 3) = vortg(iz, iy, ix, 3)

            enddo
        enddo
    enddo

    call init_fft

    call vor2vel(vortg, velog, velgradg)

    error = maxval(dabs(velog(0:nz, :, :, :) - velog_ref))

    print *, error
    stop

    call print_result_dp('Test inversion (vorticity)', error, atol=2.0e-14)

    deallocate(velog_ref)

end program test_vor2vel_2
