! =============================================================================
!                             Test tri-inversion module
!
!                     This unit test checks the inversion.
! =============================================================================
program test_tri_inversion
    use unit_test
    use constants, only : pi, twopi, one, zero, two, four, f12, f32
    use parcel_container
    use tri_inversion
    use parameters, only : extent, lower, update_parameters, vcell, dx, nx, nz, hl, ncell
    use timer
    implicit none

    double precision, allocatable :: ue(:, :), we(:, :)
    double precision, allocatable :: vortg(:, :, :), velog(:, :, :), velgradg(:, :, :)
    double precision              :: k, m, px, xx, az, mz, zz, uea, a, max_err
    integer                       :: ix, iz

    nx = 40
    nz = 20
    lower  = (/-two, -one/)
    extent = (/four, two/)

    call register_timer('vor2vel', vor2vel_timer)

    call update_parameters

    k = twopi / extent(1)
    m = pi / extent(2)

    allocate(ue(0:nz,0:nx-1))
    allocate(we(0:nz,0:nx-1))
    allocate(vortg(-1:nz+1,0:nx-1, 1))
    allocate(velog(-1:nz+1, 0:nx-1, 2))
    allocate(velgradg(-1:nz+1, 0:nx-1, 4))

    velog = zero
    vortg = zero
    velgradg = zero

    ! Set up analytical test:
    px = one
    !px = zero  !uncomment to have no mean component of the vorticity
    mz = m * lower(2)
    a = dsqrt(m**2 + k**2)
    uea = px * m * dexp(mz) * dsin(mz)     !Ensures mean ue = 0 at z = lower(2)
    do ix = 0,nx-1
        xx = k * (dx(1) * dble(ix) - hl(1))     !xx = k*x
        do iz = 0,nz
            zz = lower(2) + dx(2) * dble(iz)
            az = a * zz                  !az = a * z
            mz = m * zz                  !mz = m * z
            ue(iz,ix) = dexp(az) * (m * dsin(mz) - a * dcos(mz)) * dsin(xx) &
                      + px * dexp(mz) * m * (dsin(mz) - dcos(mz)) - uea
            we(iz,ix) = k * dexp(az) * dcos(mz) * dcos(xx)
            vortg(iz,ix, 1) = -two * m * (a * dexp(az) * dsin(xx) + px * m * dexp(mz)) * dsin(mz)
        enddo
    enddo

    call init_inversion

    call vor2vel(vortg, velog, velgradg)

    max_err = zero

    ! absolute error (reference obtained from David's program (tri_inversion.f90)
    max_err = max(max_err, abs(maxval(abs(velog(0:nz, :, 1) - ue)) - 0.72528494909253993d0))
    max_err = max(max_err, abs(maxval(abs(velog(0:nz, :, 2) - we)) - 0.000039197805283164300d0))


    ! rms error
    ue = (velog(0:nz, :, 1) - ue) ** 2
    xx = dsqrt((f12 * sum(ue(0, :) + ue(nz,:)) + sum(ue(1:nz-1, :))) / dble(ncell))
    we = (velog(0:nz, :, 2) - we) ** 2
    zz = dsqrt(sum(we(1:nz-1, :)) / dble(ncell))

    ! referene obtaind from David's program (tri_inversion.f90)
    max_err = max(max_err, abs(xx - 0.1203711977266676d0))
    max_err = max(max_err, abs(zz - 0.0000181041574751d0))

    call print_result_dp('Test tri-inversion', max_err)

    deallocate(ue)
    deallocate(we)
    deallocate(vortg)
    deallocate(velog)
    deallocate(velgradg)

end program test_tri_inversion
