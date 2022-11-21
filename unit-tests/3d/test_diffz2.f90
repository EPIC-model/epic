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
    use constants, only : zero, one, two, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils, only : init_inversion, central_diffz, diffz &
                              , field_combine_physical, field_decompose_physical
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), fp(:, :, :), &
                                     ds(:, :, :), dp(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z

    ! doubling nx, ny, nz reduces the error by a factor 16
    nx = 64
    ny = 64
    nz = 32
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
                z = lower(3) + iz * dx(3)

                fp(iz, iy, ix) = dcos(x) * dcos(y) * dsin(z)
                ref_sol(iz, iy, ix) = dcos(x) * dcos(y) * dcos(z)
            enddo
        enddo
    enddo

    call init_inversion

    call central_diffz(fp, dp)

    error = maxval(dabs(dp - ref_sol))

    print *, "central", error

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)

    error = maxval(dabs(dp - ref_sol))

    print *, "exact", error

    call print_result_dp('Test diffz', error, atol=5.0e-4)

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

end program test_diffz2
