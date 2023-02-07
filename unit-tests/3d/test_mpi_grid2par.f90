! =============================================================================
!                       Test grid2par interpolation
!
!                       This unit test checks grid2par
! =============================================================================
program test_mpi_grid2par
    use unit_test
    use mpi_communicator
    use mpi_collectives
    use mpi_layout
    use constants, only : pi, zero, one, two, three, four, five, f12, f23
    use parcel_container
    use parcel_interpl, only : grid2par, grid2par_timer
    use parcel_ellipsoid, only : get_abc
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz
    use fields, only : velog, vortg, velgradg, field_alloc
    use mpi_timer
    implicit none

    double precision              :: error
    integer                       :: ix, iy, iz, i, j, k, l, n_per_dim
    double precision              :: im, corner(3)
    double precision, allocatable :: vel(:, :), vortend(:, :), vgrad(:, :)
    logical                       :: passed = .true.

    call mpi_comm_initialise

    passed = (passed .and. (comm%err == 0))

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-one, -one, -one/)
    extent =  (/two, two, two/)

    call register_timer('grid2par', grid2par_timer)

    call update_parameters

    call field_alloc

    n_per_dim = 3

    n_parcels = n_per_dim ** 3 * nz * (box%hi(1) - box%lo(1) + 1) * (box%hi(2) - box%lo(2) + 1)
    call parcel_alloc(n_parcels)

    allocate(vel(3, n_parcels))
    allocate(vortend(3, n_parcels))
    allocate(vgrad(5, n_parcels))


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

    velog(:, :, :, 1) = one
    velog(:, :, :, 2) = two
    velog(:, :, :, 3) = three

    vortg(:, :, :, 1) = one
    vortg(:, :, :, 2) = two
    vortg(:, :, :, 3) = three

    do l = 1, 5
        velgradg(:, :, :, l) = dble(l)
    enddo

    ! we cannot check vortend since parcels%vorticity is used in grid2par
    call grid2par(vel, vortend, vgrad)

    error = zero

    do l = 1, 3
        error = max(error, maxval(dabs(vel(l, 1:n_parcels) - dble(l))))
    enddo

    do l = 1, 5
        error = max(error, maxval(dabs(vgrad(l, 1:n_parcels) - dble(l))))
    enddo

    call mpi_blocking_reduce(error, MPI_MAX)

    passed = (passed .and. (error < dble(1.0e-14)))

    deallocate(vel)
    deallocate(vortend)
    deallocate(vgrad)

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI grid2par', passed)
    endif

end program test_mpi_grid2par
