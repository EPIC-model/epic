! =============================================================================
!                       Test parcel initialisation
!
!         This unit test checks the parcel initialisation from fields.
! =============================================================================
program test_parcel_init
    use unit_test
    use constants, only : pi, zero, one, two, four, f12
    use parcel_container
    use parcel_init, only : gen_parcel_scalar_attr
    use parcel_interpl, only : par2grid
    use parcel_ellipse, only : get_ab
    use fields, only : tbuoyg, field_default
    use parameters, only : update_parameters, dx, ncell, nx, nz, lower
    use options, only : box
    implicit none

    double precision  :: xg, zg, facx, facz, argx, argz, v0
    integer:: i, ix, iz, mx, mz
    double precision :: rms, rmserr, error
    double precision, allocatable :: workg(:, :)
    double precision :: tol = 1.0d-9

     !Number of parcels per grid box = nbgx*nbgz
    integer, parameter :: nbgx = 2, nbgz = nbgx

    !Fractions of grid cell used for placing parcels:
    double precision, parameter :: dxf = one / dble(nbgx), &
                                   dzf = one / dble(nbgz)

    box%ncells = (/64, 32/)
    box%extent = (/8.0d0, four/)
    box%origin = (/-four, -two/)

    call update_parameters

    call field_default

    !Maximum number of parcels:

    !Total number of parcels:
    n_parcels = nbgx * nbgz * ncell
    call parcel_alloc(n_parcels)


    !--------------------------------------------------------
    ! Define a gridded field "tbuoyg" (this can be arbitrary):
    facx = two * pi / box%extent(1)
    facz = one / box%extent(2)
    do ix = 0, nx-1
        xg = dx(1) * dble(ix)
        argx = facx * xg
        do iz = 0, nz
            zg = dx(2) * dble(iz)
            argz = facz * zg
            tbuoyg(iz,ix) = dexp(argz) * dsin(argx+one) / (one + f12 * dcos(argx))
        enddo
    enddo

    !---------------------------------------------------------
    !Initialise parcel volume positions and volume fractions:
    v0 = dxf * dzf !use equal volume fraction for each parcel
    i = 0
    do iz = 0, nz-1
        do ix = 0, nx-1
            do mz = 1,nbgz
                do mx = 1,nbgx
                    i = i + 1
                    parcels%position(i, 1) = lower(1) + dx(1) * (dble(ix) + dxf * (dble(mx) - f12))
                    parcels%position(i, 2) = lower(2) + dx(2) * (dble(iz) + dzf * (dble(mz) - f12))
                    parcels%volume(i) = v0
                    parcels%B(i, 2) = zero
                    parcels%B(i, 1) = get_ab(v0)
                enddo
            enddo
        enddo
    enddo


    !---------------------------
    ! Generate parcel attribute:
    call gen_parcel_scalar_attr(tbuoyg, tol, parcels%buoyancy)

    !
    ! check result
    !
    error = zero

    ! Copy buoyg since it is overwritten in par2grid after:
    allocate(workg(-1:nz+1, 0:nx-1))
    workg = tbuoyg

    call par2grid

    ! Compute max and rms errors:

    ! Rms value of original field
    rms = dsqrt((f12 * sum(workg(0,:) ** 2 + workg(nz,:) ** 2) + &
                       sum(workg(1:nz-1,:) ** 2)) / dble(ncell))

    workg = tbuoyg - workg

    ! Rms error in reconstruction
    rmserr = dsqrt((f12 * sum(workg(0,:) ** 2 + workg(nz,:) ** 2) + &
                          sum(workg(1:nz-1,:) ** 2)) / dble(ncell))

    error = max(error, rmserr)

    ! Relative rms error
    error = max(error, rmserr / rms)

    ! Maximum error
    error = max(error, maxval(dabs(workg(0:nz,:))))

    call print_result_dp('Test parcel initialisation', error, atol=two * tol)

    call parcel_dealloc

    deallocate(workg)

end program test_parcel_init
