! =============================================================================
!                       Test the velocity gradient tensor
!
!  This unit test checks the calculation of the velocity gradient tensor using the
!  incompressible Taylor-Green flow:
!               u = A * cos(a * x) * sin(b * y) * sin(c * z)
!               v = B * sin(a * x) * cos(b * y) * sin(c * z)
!               w = C * sin(a * x) * sin(b * y) * cos(c * z)
!  with A = 1, a = 2, B = -1, b = 1, C = 1, c = -1; hence, the velocity
! gradient tensor is given by
!           du/dx = - a * A * sin(a * x) * sin(b * y) * sin(c * z)
!           du/dy =   b * A * cos(a * x) * cos(b * y) * sin(c * z)
!           du/dz =   c * A * cos(a * x) * sin(b * y) * cos(c * z)

!           dv/dx =   a * B * cos(a * x) * cos(b * y) * sin(c * z)
!           dv/dy = - b * B * sin(a * x) * sin(b * y) * sin(c * z)
!           dv/dz =   c * B * sin(a * x) * cos(b * y) * cos(c * z)

!           dw/dx =   a * C * cos(a * x) * sin(b * y) * cos(c * z)
!           dw/dy =   b * C * sin(a * x) * cos(b * y) * cos(c * z)
!           dw/dz = - c * C * sin(a * x) * sin(b * y) * sin(c * z)
! The vorticity is given by
!              xi = (b * C - c * B) * sin(a * x) * cos(b * y) * cos(c * z)
!             eta = (c * A - a * C) * cos(a * x) * sin(b * y) * cos(c * z)
!            zeta = (a * B - b * A) * cos(a * x) * cos(b * y) * sin(c * z)
!
!  Reference:
!  16 November 2021
!  https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex
! =============================================================================
program test_mpi_vel2vgrad
    use unit_test
    use constants, only : zero, one, two, four, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, velgradg, field_default
    use inversion_utils, only : init_inversion
    use inversion_mod, only : vel2vgrad
    use sta3dfft, only : fftxyp2s, finalise_fft
    use mpi_timer
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: strain(:, :, :, :), svelog(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, AA, BB, CC, a, b, c
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 32
    ny = 64
    nz = 128
    lower  = (/-pi, -twopi, -two * twopi/)
    extent =  (/twopi, two * twopi, four * twopi/)

    AA =  one
    a  =  two
    BB = -one
    b  =  one
    CC =  one
    c  = -one

    call update_parameters

    call field_default

    allocate(strain(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 5))
    allocate(svelog(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

    call init_inversion

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                ! velocity
                velog(iz, iy, ix, 1) = AA * dcos(a * x) * dsin(b * y) * dsin(c * z)
                velog(iz, iy, ix, 2) = BB * dsin(a * x) * dcos(b * y) * dsin(c * z)
                velog(iz, iy, ix, 3) = CC * dsin(a * x) * dsin(b * y) * dcos(c * z)

                ! vorticity
                vortg(iz, iy, ix, 1) = (b * CC - c * BB) * dsin(a * x) * dcos(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 2) = (c * AA - a * CC) * dcos(a * x) * dsin(b * y) * dcos(c * z)
                vortg(iz, iy, ix, 3) = (a * BB - b * AA) * dcos(a * x) * dcos(b * y) * dsin(c * z)

                ! velocity gradient tensor (reference solution)
                strain(iz, iy, ix, 1) = - a * AA * dsin(a * x) * dsin(b * y) * dsin(c * z) ! du/dx
                strain(iz, iy, ix, 2) =   b * AA * dcos(a * x) * dcos(b * y) * dsin(c * z) ! du/dy
                strain(iz, iy, ix, 3) = - b * BB * dsin(a * x) * dsin(b * y) * dsin(c * z) ! dv/dy
                strain(iz, iy, ix, 4) =   a * CC * dcos(a * x) * dsin(b * y) * dcos(c * z) ! dw/dx
                strain(iz, iy, ix, 5) =   b * CC * dsin(a * x) * dcos(b * y) * dcos(c * z) ! dw/dy
            enddo
        enddo
    enddo

    call fftxyp2s(velog(:, :, :, 1), svelog(:, :, :, 1))
    call fftxyp2s(velog(:, :, :, 2), svelog(:, :, :, 2))
    call fftxyp2s(velog(:, :, :, 3), svelog(:, :, :, 3))

    call vel2vgrad(svelog)

    call finalise_fft

    error = maxval(dabs(velgradg(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :)   &
                        - strain(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :)))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(error, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%master, comm%world, comm%err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0) .and. (error < dble(2.0e-14)))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI vel2vgrad', passed)
    endif

    deallocate(strain)
    deallocate(svelog)

end program test_mpi_vel2vgrad
