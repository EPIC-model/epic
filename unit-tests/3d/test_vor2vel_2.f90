! =============================================================================
!                       Test the vorticity inversion !FIXME This test is now the TG flow
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
    use constants, only : zero, one, two, four, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, velgradg, field_alloc
    use inversion_utils, only : init_fft, fftxyp2s
    use inversion_mod, only : vor2vel, vor2vel_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: velog_ref(:, :, :, :), ozm(:, :), svelog(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, AA, BB, CC, a, b, c

    call register_timer('vorticity', vor2vel_timer)


    nx = 32
    ny = 32
    nz = 32
    lower  = (/-pi, -pi, -pi/)
    extent =  (/twopi, twopi, twopi/)

    AA =  one
    a  =  two
    BB = -one
    b  =  one
    CC =  one
    c  = -one

    allocate(velog_ref(0:nz, 0:ny-1, 0:nx-1, 3))
    allocate(ozm(0:ny-1, 0:nx-1))
!     allocate(svelog(0:nz, 0:nx-1, 0:ny-1, 3))

    call update_parameters

    call field_alloc

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                ! vorticity
                vortg(iz, iy, ix, 1) = (b * CC - c * BB) * dsin(a * x) * dcos(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 2) = (c * AA - a * CC) * dcos(a * x) * dsin(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 3) = (a * BB - b * AA) * dcos(a * x) * dcos(b * y) * dsin(c * z)

                ! velocity (reference solution)
                velog_ref(iz, iy, ix, 1) = AA * dcos(a * x) * dsin(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 2) = BB * dsin(a * x) * dcos(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 3) = CC * dsin(a * x) * dsin(b * y) * dcos(c * z)

!                 print *, x, y, z, vortg(iz, iy, ix, 1), vortg(iz, iy, ix, 2), vortg(iz, iy, ix, 3)
!                 print *, x, y, z, velog_ref(iz, iy, ix, 1), velog_ref(iz, iy, ix, 2), velog_ref(iz, iy, ix, 3)
            enddo
        enddo
    enddo

!     stop

    call init_fft

!     ozm = (f12 * (vortg(-1, :, :, 3) + vortg(nz+1, :, :, 3)) &
!         + sum(vortg(1:nz, :, :, 3))) / dble(nz+1)
!
!     do iz = -1, nz+1
!         vortg(iz, :, :, 3) = vortg(iz, :, :, 3) - ozm
!     enddo

    call vor2vel(vortg, velog, velgradg)

!     velog(0:nz, :, :, :) = dabs(velog(0:nz, :, :, :) - velog_ref)
!
!     do ix = 0, nx-1
!         x = lower(1) + ix * dx(1)
!         do iy = 0, ny-1
!             y = lower(2) + iy * dx(2)
!             do iz = 0, nz
!                 z = lower(3) + iz * dx(3)
!                 print *, x, y, z, velog(iz, iy, ix, 1), velog(iz, iy, ix, 2), velog(iz, iy, ix, 3)
!             enddo
!         enddo
!     enddo
!     stop


    error = maxval(dabs(velog(0:nz, :, :, :) - velog_ref))

    print *, error

    call print_result_dp('Test inversion (vorticity)', error, atol=2.0e-14)

    deallocate(velog_ref)
    deallocate(ozm)

end program test_vor2vel_2
