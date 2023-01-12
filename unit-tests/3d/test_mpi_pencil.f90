! =============================================================================
!                       Test pencil module
!
!       This unit test checks the pencil transposition.
! =============================================================================
program test_mpi_pencil
    use unit_test
    use constants, only : zero, one, two, four
    use sta2dfft
    use pencil_fft
    use dimensions
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    implicit none

    double precision, allocatable :: values(:, :, :), vtrans(:, :, :)
    logical                       :: passed
    integer                       :: iy !ix, iy, iz

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    nx = 8
    ny = 16
    nz = 32
    lower = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(vtrans(box%hlo(2):box%hhi(2), box%hlo(3):box%hhi(3), box%hlo(1):box%hhi(1)))

    values(:, :, :) = zero
    vtrans(:, :, :) = zero

!     print *, box%size
!     stop

!     do iz = box%lo(3), box%hi(3)
!         do ix = box%lo(1), box%hi(1)
            do iy = box%lo(2), box%hi(2)
                values(:, iy, :) = iy+1
            enddo
!             write(WRITE_VOR, '(1x,f13.6,6(1x,1p,e14.7))')
!             write(*, '(8(f1.0))')
!             print *, values(iz, 0:ny-1, ix)
!         enddo
!         print *, ""
!     enddo

    call initialise_pencil_fft(nx, ny, nz)


    call transpose_to_pencil(y_from_z_transposition, (/1, 2, 3/), dim_y_comm, FORWARD,        &
                             values(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                             fft_in_y_buffer)

    if (mpi_rank == mpi_master) then
!         do ix = box%lo(1), box%hi(1)
!             do iy = box%lo(2), box%hi(2)
                do iy = 1, ny
                    print *, fft_in_y_buffer(iy, 1, 1)
                enddo
!             enddo
!         enddo
    endif

    call finalise_pencil_fft


    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test pencil transposition', passed)
    endif


end program test_mpi_pencil
