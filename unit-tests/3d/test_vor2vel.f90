! =============================================================================
!                       Test the vorticity inversion
!
!  This unit test checks the calculation of the velocity using the
!  Taylor-Green flow:
!               u = A * cos(a * x) * sin(b * y) * sin(c * z)
!               v = B * sin(a * x) * cos(b * y) * sin(c * z)
!               w = C * sin(a * x) * sin(b * y) * cos(c * z)
!  The vorticity of this flow is
!               xi = (b * C - c * B) * sin(a * x) * cos(b * y) * cos(c * z)
!              eta = (c * A - a * C) * cos(a * x) * sin(b * y) * cos(c * z)
!             zeta = (a * B - b * A) * cos(a * x) * cos(b * y) * sin(c * z)
! =============================================================================
program test_vor2vel_2
    use unit_test
    use constants, only : zero, one, two, four, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, field_alloc
    use inversion_utils, only : init_inversion
    use inversion_mod, only : vor2vel, vor2vel_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: velog_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, AA, BB, CC, a, b, c

    call register_timer('vorticity', vor2vel_timer)

    nx = 64
    ny = 64
    nz = 64

    lower  = (/-pi, -pi, -f12 * pi/)
    extent =  (/twopi, twopi, twopi/)

    AA = one
    BB = one
    CC = -two
    a = one
    b = one
    c = one

    allocate(velog_ref(-1:nz+1, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_alloc

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = -1, nz+1
                z = lower(3) + iz * dx(3)

                ! vorticity
                vortg(iz, iy, ix, 1) = (b * CC - c * BB) * dsin(a * x) * dcos(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 2) = (c * AA - a * CC) * dcos(a * x) * dsin(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 3) = (a * BB - b * AA) * dcos(a * x) * dcos(b * y) * dsin(c * z)

                ! velocity (reference solution)
                velog_ref(iz, iy, ix, 1) = AA * dcos(a * x) * dsin(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 2) = BB * dsin(a * x) * dcos(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 3) = CC * dsin(a * x) * dsin(b * y) * dcos(c * z)

            enddo
        enddo
    enddo

    call init_inversion

    call vor2vel

    error = maxval(dabs(velog_ref(0:nz, :, :, :) - velog(0:nz, :, :, :)))

    call print_result_dp('Test inversion (vorticity)', error, atol=1.0e-14)

    deallocate(velog_ref)

end program test_vor2vel_2
