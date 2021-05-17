! =============================================================================
!                       Test trilinear interpolation
!
!         This unit test checks the trilinear interoplation par2grid
! =============================================================================
program test_trilinear
    use unit_test
    use constants, only : pi, zero, one
    use parcel_container
    use parcel_interpl, only : par2grid, grid2par
    use options, only : parcel_info, box, interpl
    use parameters, only : lower, update_parameters, vcell, dx, nx, nz, ngrid
    implicit none

    double precision :: volg(-1:33, 0:31), error
    integer :: i, j, k, jj, ii

    box%nc = (/32, 32/)
    box%extent =  (/0.4d0, 0.4d0/)

    call update_parameters()

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)

    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + 0.25d0 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + 0.25d0 * dx(2) * jj
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    volg = zero

    parcel_info%is_elliptic = .true.

    parcels%volume = 0.25d0 * vcell

    ! b11
    parcels%B(:, 1) = 0.25d0 * vcell / pi

    ! b12
    parcels%B(:, 2) = zero


    interpl = 'trilinear'

    call par2grid(parcels)

    error = abs(sum(volg(0:nz, 0:nx-1)) - dble(ngrid) * vcell)

    call print_result_dp('Test trilinear (par2grid)', error, atol=dble(1.0e-14))

end program test_trilinear
