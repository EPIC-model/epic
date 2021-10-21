! =============================================================================
!                       Test linterp interpolation
!
!         This unit test checks the linterp interpolation par2grid
! =============================================================================
program test_linterp
    use unit_test
    use constants, only : pi, zero, one, f14, f32
    use parcel_container
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipse, only : get_ab
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : volg, field_alloc
    use timer
    implicit none

    double precision :: error
    integer :: i, j, k, l, kk, jj, ii

    ny = 32
    nz = 32
#ifdef ENABLE_3D
    nx = 32
    lower  = (/-f32, -f32, -f32/)
    extent =  (/0.4d0, 0.4d0, 0.4d0/)
#else
    nx = 0
    lower  = (/-f32, -f32/)
    extent =  (/0.4d0, 0.4d0/)
#endif

    call register_timer('par2grid', par2grid_timer)

    call update_parameters

    call field_alloc

    n_parcels = 4*nx*ny*nz
    call parcel_alloc(n_parcels)

    l = 1
    do k = 0, nz-1
        do j = 0, ny-1
            do i = 0, nx-1
                do kk = 1, 4, 2
                    do jj = 1, 4, 2
                        do ii = 1, 4, 2
                            parcels%position(l, 1) = lower(1) + i * dx(1) + f14 * dx(1) * ii
                            parcels%position(l, 2) = lower(2) + j * dx(2) + f14 * dx(2) * jj
#ifdef ENABLE_3D
                            parcels%position(l, 3) = lower(3) + k * dx(3) + f14 * dx(3) * kk
#endif
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    parcels%volume = f14 * vcell

    ! b11
    parcels%B(:, 1) = get_ab(parcels%volume(1:n_parcels))

    ! b12
    parcels%B(:, 2) = zero

    parcels%B(:, 3:5) = zero

    call par2grid

    error = abs(sum(volg(0:nz, 0:ny-1, 0:nx-1)) - dble(ngrid) * vcell)

    call print_result_dp('Test linterp (par2grid)', error, atol=dble(1.0e-14))

end program test_linterp
