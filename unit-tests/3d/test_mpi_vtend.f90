! =============================================================================
!                       Test the vorticity tendency
!
!  This unit test checks the calculation of the vorticity tendency using the
!  Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
!  The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
!  For this setting the theoretical vorticity tendency is
!  vtend(z, y, x, 1) = alpha*k*m**2 * (k^2+l^2)^(-1)*sin(k*x+l*y)*cos(k*x+l*y)
!  vtend(z, y, x, 2) = alpha*l*m**2 * (k^2+l^2)^(-1)*sin(k*x+l*y)*cos(k*x+l*y)
!  vtend(z, y, x, 3) = -alpha * m * sin(m*z) * cos(m*z)
! =============================================================================
program test_mpi_vtend
    use unit_test
    use constants, only : one, two, pi, f12, f34, three
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields, only : vortg, velog, vtend, field_default
    use sta3dfft, only : finalise_fft, fftxyp2s
    use inversion_mod, only : init_inversion, vorticity_tendency, vtend_timer
    use timer
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: vtend_ref(:, :, :, :)
    integer                       :: ix, iy, iz, ik, il, im
    double precision              :: x, y, z, alpha, fk2l2, k, l, m
    double precision              :: cosmz, sinmz, sinkxly, coskxly
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    call register_timer('vtend', vtend_timer)

    nx = 128
    ny = 128
    nz = 128

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call update_parameters

    call field_default

    allocate(vtend_ref(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 3))

    call init_inversion

    do ik = 1, 3
        k = dble(ik)
        do il = 1, 3
            l = dble(il)
            do im = 1, 3, 2
                m = dble(im)

                alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
                fk2l2 = one / dble(k ** 2 + l ** 2)

                do ix = box%lo(1), box%hi(1)
                    x = lower(1) + ix * dx(1)
                    do iy = box%lo(2), box%hi(2)
                        y = lower(2) + iy * dx(2)
                        do iz = -1, nz+1
                            z = lower(3) + iz * dx(3)

                            cosmz = dcos(m * z)
                            sinmz = dsin(m * z)
                            sinkxly = dsin(k * x + l * y)
                            coskxly = dcos(k * x + l * y)

                            ! velocity
                            velog(iz, iy, ix, 1) = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                            velog(iz, iy, ix, 2) = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                            velog(iz, iy, ix, 3) = cosmz * coskxly

                            ! vorticity
                            vortg(iz, iy, ix, 1) = alpha * velog(iz, iy, ix, 1)
                            vortg(iz, iy, ix, 2) = alpha * velog(iz, iy, ix, 2)
                            vortg(iz, iy, ix, 3) = alpha * velog(iz, iy, ix, 3)

                            ! reference solution
                            vtend_ref(iz, iy, ix, 1) = alpha * k * m ** 2 * fk2l2 * sinkxly * coskxly
                            vtend_ref(iz, iy, ix, 2) = alpha * l * m ** 2 * fk2l2 * sinkxly * coskxly
                            vtend_ref(iz, iy, ix, 3) = -alpha * m * sinmz * cosmz
                        enddo
                    enddo
                enddo

                call vorticity_tendency

                error = max(error, maxval(dabs(vtend_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :) &
                                                 - vtend(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), :))))
            enddo
        enddo
    enddo

    call finalise_fft

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

    passed = (passed .and. (comm%err == 0) .and. (error < 1.2d0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test vorticity tendency', passed)
    endif

    deallocate(vtend_ref)

end program test_mpi_vtend
