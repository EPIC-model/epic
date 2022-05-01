! =============================================================================
!                       Test parcel initialisation in 3D
!
!         This unit test checks the parcel initialisation from fields.
! =============================================================================
program test_mpi_parcel_init_3d
    use unit_test
    use mpi_communicator
    use field_mpi
    use constants, only : pi, zero, one, two, four, five, f12, f13, f23, f32
    use parcel_container
    use parcel_init, only : gen_parcel_scalar_attr, unit_test_parcel_init_alloc, init_timer
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipsoid, only : get_abc
    use fields, only : tbuoyg, field_default
    use field_ops, only : get_rms, get_abs_max
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, vcell
    use timer
    implicit none

    double precision  :: xg, yg, zg, facx, facy, facz, argx, argy, argz, v0
    integer :: i, j, k, ix, iy, iz, l
    double precision :: rms, rmserr, error, corner(3), im(3)
    logical          :: passed = .true.
    double precision, allocatable :: workg(:, :, :)
    double precision :: tol = 1.0d-9

     !Number of parcels per grid box = nbgx*nbgz
    integer, parameter :: nbgx = 2, nbgy = nbgx, nbgz = nbgx

    !Fractions of grid cell used for placing parcels:
    double precision, parameter :: dxf = one / dble(nbgx), &
                                   dyf = one / dble(nbgy), &
                                   dzf = one / dble(nbgz)

    call mpi_comm_initialise

    call parse_command_line

    passed = (mpi_err == 0)

    nx = 64
    ny = 64
    nz = 32
    lower = (/-four, -four, -two/)
    extent = (/8.0d0, 8.0d0, four/)

    call register_timer('parcel init', init_timer)
    call register_timer('par2grid', par2grid_timer)

    call update_parameters

    call field_default

    !Maximum number of parcels:

    !Total number of parcels:
    n_parcels = nbgx * nbgy * nbgz * nz * (box%hi(1) - box%lo(1) + 1) * (box%hi(2) - box%lo(2) + 1)
    call parcel_alloc(n_parcels)

    !--------------------------------------------------------
    ! Define a gridded field "tbuoyg" (this can be arbitrary):
    facx = two * pi / extent(1)
    facy = two * pi / extent(2)
    facz = one / extent(3)
    do ix = box%lo(1), box%hi(1)
        xg = dx(1) * dble(ix)
        argx = facx * xg
        do iy = box%lo(2), box%hi(2)
            yg = dx(2) * dble(iy)
            argy = facy * yg
            do iz = 0, nz
                zg = dx(3) * dble(iz)
                argz = facz * zg
                tbuoyg(iz, iy, ix) = dexp(argz) * dsin(argx - argy + one) / &
                                     ((one + f12 * dcos(argx)) * (one + f13 * dcos(argy)))
            enddo
        enddo
    enddo

    call field_halo_swap(tbuoyg)

    !---------------------------------------------------------
    !Initialise parcel volume positions and volume fractions:
    v0 = dxf * dyf * dzf * vcell
    l = 1

    im(1) = one / dble(nbgx)
    im(2) = one / dble(nbgy)
    im(3) = one / dble(nbgz)

    do iz = 0, nz-1
        do iy = box%lo(2), box%hi(2)
            do ix = box%lo(1), box%hi(1)
                corner = lower + dble((/ix, iy, iz/)) * dx
                do k = 1, nbgz
                    do j = 1, nbgy
                        do i = 1, nbgx
                            parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im(1)
                            parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im(2)
                            parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im(3)
                            parcels%volume(l) = v0
                            parcels%B(:, l) = zero
                            parcels%B(1, l) = get_abc(v0) ** f23
                            parcels%B(4, l) = parcels%B(1, l)
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    ! Prepare for "gen_parcel_scalar_attr"
    call unit_test_parcel_init_alloc

    !---------------------------
    ! Generate parcel attribute:
    call gen_parcel_scalar_attr(tbuoyg, tol, parcels%buoyancy)

    !
    ! check result
    !
    error = zero

    ! Copy buoyg since it is overwritten in par2grid after:
    allocate(workg(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    workg = tbuoyg

    call par2grid

    ! Compute max and rms errors:

    ! Rms value of original field
    rms = get_rms(workg)

    workg = tbuoyg - workg

    ! Rms error in reconstruction
    rmserr = get_rms(workg)

    error = max(error, rmserr)

    ! Relative rms error
    error = max(error, rmserr / rms)

    ! Maximum error
    error = max(error, get_abs_max(workg))

    passed = (passed .and. (error < two * tol))

    call parcel_dealloc

    deallocate(workg)

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test parcel initialisation 3D', passed)
    endif

end program test_mpi_parcel_init_3d
