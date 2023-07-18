! =============================================================================
!                       Test trilinear interpolation
!
!         This unit test checks the trilinear interpolation par2grid
! =============================================================================
program test_mpi_trilinear
    use unit_test
    use mpi_communicator
    use constants, only : pi, zero, one, f12, f23, f32
    use parcel_container
    use mpi_layout
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipsoid, only : get_abc
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : volg, field_alloc
    use field_ops, only : get_sum
    use mpi_timer
    implicit none

    double precision :: error
    integer          :: ix, iy, iz, i, j, k, l, n_per_dim
    double precision :: im, corner(3)
    logical          :: passed = .true.

    call mpi_comm_initialise

    passed = (passed .and. (world%err == 0))

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-f32, -f32, -f32/)
    extent =  (/0.4d0, 0.4d0, 0.4d0/)

    call register_timer('par2grid', par2grid_timer)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_alloc

    n_per_dim = 3

    n_parcels = n_per_dim ** 3 * nz * (box%hi(1) - box%lo(1) + 1) * (box%hi(2) - box%lo(2) + 1)
    call parcel_alloc(n_parcels)


    im = one / dble(n_per_dim)

    l = 1
    do iz = 0, nz-1
        do iy = box%lo(2), box%hi(2)
            do ix = box%lo(1), box%hi(1)
                corner = lower + dble((/ix, iy, iz/)) * dx
                do k = 1, n_per_dim
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                            parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    parcels%volume = vcell / dble(n_per_dim ** 3)

    parcels%B(:, 1:n_parcels) = zero

    ! b11
    parcels%B(1, 1:n_parcels) = get_abc(parcels%volume(1:n_parcels)) ** f23

    ! b22
    parcels%B(4, 1:n_parcels) = parcels%B(1, 1:n_parcels)

    call par2grid

    error = abs(get_sum(volg) - dble(ngrid) * vcell)

    passed = (passed .and. (error < dble(3.0e-14)))

    call mpi_comm_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI trilinear (par2grid)', passed)
    endif

end program test_mpi_trilinear
