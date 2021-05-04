! =============================================================================
!                       Test trilinear interpolation
!
!         This unit test checks the trilinear interoplation par2grid
! =============================================================================
program test_trilinear
    use constants, only : pi
    use parcel_container
    use parcel_interpl, only : par2grid, grid2par
    use options, only : parcel_info, grid, interpl
    use parameters, only : extent, lower, update_parameters, vcell, dx, nx, nz, ngrid
    implicit none

    double precision :: volg(-1:33, 0:31, 1), error
    integer :: i, j, k, jj, ii

    grid = (/33, 33/)

    extent =  (/0.4, 0.4/)
    lower = (/-0.2, -0.2/)

    call update_parameters()

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)


    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + 0.25 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + 0.25 * dx(2) * jj
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    volg = 0.0

    parcel_info%is_elliptic = .true.

    parcels%volume = 0.25d0 * vcell

    ! b11
    parcels%B(:, 1) = 1.0d0

    ! b12
    parcels%B(:, 2) = 0.0d0


    interpl = 'trilinear'

    call par2grid(parcels, parcels%volume, volg)

    error = abs(sum(volg(0:nz, 0:nx-1, :)) - ngrid * vcell)

    if (error > 1.0e-15) then
        print '(a27, a9)', 'Test trilinear (par2grid):', 'FAILED'
    else
        print '(a27, a9)', 'Test trilinear (par2grid):', 'PASSED'
    endif

end program test_trilinear
