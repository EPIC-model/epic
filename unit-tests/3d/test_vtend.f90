! =============================================================================
!                       Test the vorticity tendency
!
!  This unit test checks the calculation of the vorticity tendency using the
!  ABC flow:
!               u = A * sin(z) + C * cos(y)
!               v = B * sin(x) + A * cos(z)
!               w = C * sin(y) + B * cos(x)
!  The vorticity of this flow is
!               xi = A * sin(z) + C * cos(y)
!              eta = B * sin(x) + A * cos(z)
!             zeta = C * sin(y) + B * cos(x)
!  hence, the cross-product of vorticity and velocity (i.e. Lamb vector)
!  is zero. The vorticity tendency is given by
!               Sx = A * B * cos(x) * cos(z) - B * C * sin(x) * sin(y) + db/dy
!               Sy = B * C * cos(x) * cos(y) - A * C * sin(y) * sin(z) - db/dx
!               Sz = A * C * cos(y) * cos(z) - A * B * sin(x) * sin(z)
!  with buoyancy derivatives db/dy and db/dx. In this test we use a buoyancy
!  of the form:
!               b(x, y, z) = b(x, y) = cos(x) + sin(y)
!                    db/dx = - sin(x)
!                    db/dy = cos(y)
!
!  Reference:
!  16 November 2021
!  https://en.wikipedia.org/wiki/Arnold%E2%80%93Beltrami%E2%80%93Childress_flow
! =============================================================================
program test_vtend
    use unit_test
    use constants, only : zero, two, four, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velgradg, vtend, tbuoyg, field_alloc
    use inversion_utils, only : init_fft
    use inversion_mod, only : vorticity_tendency
    use timer
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: S(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, A, B, C


    nx = 32
    ny = 64
    nz = 128
    lower  = (/-pi, -twopi, -two * twopi/)
    extent =  (/twopi, two * twopi, four * twopi/)

    A = 0.5d0
    B = 1.2d0
    C = 0.9d0

    allocate(S(0:nz, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_alloc

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                ! buoyancy
                tbuoyg(iz, iy, ix) = dcos(x) + dsin(y)

                ! vorticity
                vortg(iz, iy, ix, 1) = A * dsin(z) + C * dcos(y)
                vortg(iz, iy, ix, 2) = B * dsin(x) + A * dcos(z)
                vortg(iz, iy, ix, 3) = C * dsin(y) + B * dcos(x)

                ! velocity gradient tensor
                velgradg(iz, iy, ix, 1) = zero          ! du/dx
                velgradg(iz, iy, ix, 2) = - C * dsin(y) ! du/dy
                velgradg(iz, iy, ix, 3) =   A * dcos(z) ! du/dz

                velgradg(iz, iy, ix, 4) =   B * dcos(x) ! dv/dx
                velgradg(iz, iy, ix, 5) = zero          ! dv/dy
                velgradg(iz, iy, ix, 6) = - A * dsin(z) ! dv/dz

                velgradg(iz, iy, ix, 7) = - B * dsin(x) ! dw/dx
                velgradg(iz, iy, ix, 8) =   C * dcos(y) ! dw/dy
                velgradg(iz, iy, ix, 9) = zero          ! dw/dz

                ! vorticity tendency (reference solution)
                S(iz, iy, ix, 1) = A * B * dcos(x) * dcos(z) - B * C * dsin(x) * dsin(y) + dcos(y)
                S(iz, iy, ix, 2) = B * C * dcos(x) * dcos(y) - A * C * dsin(y) * dsin(z) + dsin(x)
                S(iz, iy, ix, 3) = A * C * dcos(y) * dcos(z) - A * B * dsin(x) * dsin(z)

            enddo
        enddo
    enddo

    call init_fft

    call vorticity_tendency(vortg, tbuoyg, velgradg, vtend)

    error = maxval(dabs(vtend(0:nz, :, :, :) - S))

    print *, error

    call print_result_dp('Test inversion (vorticity tendency)', error)

    deallocate(S)

end program test_vtend
