! =============================================================================
!                             Test tri-inversion module
!
!                     This unit test checks the inversion.
! =============================================================================
program test_tri_inversion
    use unit_test
    use constants, only : pi, twopi, one, zero
    use parcel_container
    use tri_inversion
    use options, only : box
    use parameters, only : extent, lower, update_parameters, vcell, dx, nx, nz, hl, ncell
    implicit none

    double precision :: ue(0:nz,0:nx-1),we(0:nz,0:nx-1)
    double precision :: k, m, px, xx, az, mz, zz, uea, a, max_err
    integer          :: ix

    box%nc = (/40, 20/)
    box%extent = (/4.0d0, 2.0d0/)

    call update_parameters()

    k = twopi / extent(1)
    m = pi / extent(2)

    ! Set up analytical test:
    px = one
    !px = zero  !uncomment to have no mean component of the vorticity
    mz = m * zmin
    a = dsqrt(m**2 + k**2)
    uea = px * m * dexp(mz) * dsin(mz)     !Ensures mean ue = 0 at z = zmin
    do ix = 0,nx-1
        xx = k * (dx(1) * dble(ix) - hl(1))     !xx = k*x
        do iz = 0,nz
            zz =zmin + dx(2) * dble(iz)
            az = a * zz                  !az = a * z
            mz = m * zz                  !mz = m * z
            ue(iz,ix) = dexp(az) * (m * dsin(mz) - a * dcos(mz)) * dsin(xx) &
                      + px * dexp(mz) * m * (dsin(mz) - dcos(mz)) - uea
            we(iz,ix) = k * dexp(az) * dcos(mz) * dcos(xx)
            pp(iz,ix) = -two * m * (a * dexp(az) * dsin(xx) + px * m * dexp(mz)) * dsin(mz)
        enddo
    enddo

    init_error = abs(sum(volg(0:nz, 0:nx-1, :)) - ngrid * vcell)

    call init_inversion

    call vor2vel(pp, uu, ww)

    max_err = zero

    ! absolute error
    max_err = max(max_err, maxval(abs(uu(0:nz, :) - ue))
    max_err = max(max_err, maxval(abs(ww(0:nz,:) - we))

    ! rms error
    ue = (uu(0:nz, :) - ue) ** 2
    xx = dsqrt((f12 * sum(ue(0, :) + ue(nz,:)) + sum(ue(1:nz-1, :))) / dble(ncell))
    we = (ww(0:nz, :) - we) ** 2
    zz = dsqrt(sum(we(1:nz-1, :)) / dble(ncell))

    max_err = max(max_err, xx)
    max_err = max(max_err, zz)

    call print_result_dp('Test tri-inversion', max_err)

end program test_tri_inversion
