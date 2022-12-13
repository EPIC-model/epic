! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz and central_diffz using the
!  function:
!               cos(x) cos(y) cos(z)
!  in a domain of width 2 * pi in x and y, and of height pi (0 < z < pi)
!  The subroutine diffz and central_diffz should return
!               - cos(x) cos(y) sin(z)
! =============================================================================
program test_diffz2
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils, only : init_inversion, central_diffz, diffz &
                              , field_combine_physical, field_decompose_physical, fftxyp2s, fftxys2p &
                              , field_decompose_semi_spectral, field_combine_semi_spectral
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), fp(:, :, :), &
                                     ds(:, :, :), dp(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z

    nx = 128
    ny = 128
    nz = 64
    lower  = (/zero, zero, zero/)
    extent =  (/twopi, twopi, pi/)

    allocate(fp(0:nz, ny, nx))
    allocate(dp(0:nz, ny, nx))
    allocate(fs(0:nz, nx, ny))
    allocate(ds(0:nz, nx, ny))
    allocate(ref_sol(0:nz, ny, nx))

    call update_parameters

    do ix = 1, nx
        x = lower(1) + (ix - 1) * dx(1)
        do iy = 1, ny
            y = lower(2) + (iy - 1) * dx(2)
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                fp(iz, iy, ix) = dcos(x) * dcos(y) * dcos(z)
                ref_sol(iz, iy, ix) = -dcos(x) * dcos(y) * dsin(z)
            enddo
        enddo
    enddo

    call init_inversion

    dp = fp
    call fftxyp2s(dp, fs)
    call central_diffz(fs, ds)
    call fftxys2p(ds, dp)

    error = maxval(dabs(dp - ref_sol))

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)


    error = max(error, maxval(dabs(dp - ref_sol)))

    call print_result_dp('Test diffz', error, atol=4.0e-2)

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

    contains

end program test_diffz2
