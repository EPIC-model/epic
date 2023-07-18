! =============================================================================
!                       Test the vorticity inversion
!
!  This unit test checks the calculation of the velocity using the
!  Taylor-Green flow:
!               u = A * cos(a * x) * sin(b * y) * sin(c * z)
!               v = B * sin(a * x) * cos(b * y) * sin(c * z)
!               w = C * sin(a * x) * sin(b * y) * cos(c * z)
!  The vorticity of this flow is
!               xi = (b * C - c * B) * sin(a * x) * cos(b * y) * cos(c * z)
!              eta = (c * A - a * C) * cos(a * x) * sin(b * y) * cos(c * z)
!             zeta = (a * B - b * A) * cos(a * x) * cos(b * y) * sin(c * z)
! =============================================================================
program test_mpi_vor2vel
    use unit_test
    use constants, only : zero, one, two, four, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, field_default
    use inversion_utils, only : init_inversion
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: velog_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, AA, BB, CC, a, b, c
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (world%err == 0)

    call register_timer('vorticity', vor2vel_timer)

    nx = 64
    ny = 64
    nz = 64

    lower  = (/-pi, -pi, -f12 * pi/)
    extent =  (/twopi, twopi, twopi/)

    AA = one
    BB = one
    CC = -two
    a = one
    b = one
    c = one

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    allocate(velog_ref(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 3))

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = -1, nz+1
                z = lower(3) + iz * dx(3)

                ! vorticity
                vortg(iz, iy, ix, 1) = (b * CC - c * BB) * dsin(a * x) * dcos(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 2) = (c * AA - a * CC) * dcos(a * x) * dsin(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 3) = (a * BB - b * AA) * dcos(a * x) * dcos(b * y) * dsin(c * z)

                ! velocity (reference solution)
                velog_ref(iz, iy, ix, 1) = AA * dcos(a * x) * dsin(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 2) = BB * dsin(a * x) * dcos(b * y) * dsin(c * z)
                velog_ref(iz, iy, ix, 3) = CC * dsin(a * x) * dsin(b * y) * dcos(c * z)

            enddo
        enddo
    enddo

    call init_inversion

    call vor2vel

    error = maxval(dabs(velog_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :) &
                          - velog(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :)))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, world%root, world%comm, world%err)
    else
        call MPI_Reduce(error, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, world%root, world%comm, world%err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (world%err == 0) .and. (error < dble(1.0e-14)))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI vor2vel', passed)
    endif

    deallocate(velog_ref)

end program test_mpi_vor2vel
