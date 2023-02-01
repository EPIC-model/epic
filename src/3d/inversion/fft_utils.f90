module fft_utils
    use mpi_communicator
    use mpi_layout
    use constants, only : zero
    use field_mpi, only : field_halo_fill
    use parameters, only : nx, ny
    use inversion_utils, only : hrkx
    use dimensions, only : I_X, I_Y, I_Z
    implicit none

    type :: reorder_type
        integer                       :: dir                            ! reorder index
        integer, allocatable          :: dest(:)                        ! destination rank
        integer, allocatable          :: send_recv_count(:)
        integer, allocatable          :: recv_count(:)                  ! number of receives of
                                                                        ! the reordering dimension
        integer, allocatable          :: send_offset(:), recv_offset(:)
        double precision, allocatable :: send_buffer(:), recv_buffer(:)
    end type reorder_type

    logical :: l_initialised = .false.

    private :: copy_from_buffer_in_x,   &
               copy_to_buffer_in_x,     &
               copy_from_buffer_in_y,   &
               copy_to_buffer_in_y,     &
               l_initialised,           &
               setup_reordering

    type(sub_communicator) :: x_comm
    type(sub_communicator) :: y_comm

    type(reorder_type) :: x_reo
    type(reorder_type) :: y_reo

    contains

        subroutine initialise_diff

            if (l_initialised) then
                return
            endif
            l_initialised = .true.

            call setup_reordering(x_reo, x_comm, I_X, nx)

!             call setup_reordering(y_reo, y_comm, I_Y, ny)

        end subroutine initialise_diff

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine setup_reordering(reo, sub_comm, dir, ncell)
            type(reorder_type),                intent(inout) :: reo
            type(sub_communicator), intent(inout) :: sub_comm
            integer,                           intent(in)    :: dir
            integer,                           intent(in)    :: ncell
            integer, allocatable                             :: rlos(:), rhis(:)
            integer                                          :: i, j, lo, hi, rank, d

            reo%dir = dir

            ! if dir = 1 --> d = 2, i.e. if dir = I_X --> d = I_Y
            ! if dir = 2 --> d = 1, i.e. if dir = I_Y --> d = I_X
            d = mod(dir, 2) + 1

            call MPI_Cart_sub(comm%cart, (/.true., .false./), sub_comm%comm, comm%err)

            call MPI_Comm_rank(sub_comm%comm, sub_comm%rank, sub_comm%err)

            call MPI_Comm_size(sub_comm%comm, sub_comm%size, sub_comm%err)

            allocate(sub_comm%coord(1))
            call MPI_Cart_coords(sub_comm%comm, sub_comm%rank, 1, &
                                 sub_comm%coord, sub_comm%err)

            print *, comm%rank, sub_comm%rank, sub_comm%size, sub_comm%coord

            allocate(reo%send_recv_count(sub_comm%size))
            allocate(reo%send_offset(sub_comm%size))
            allocate(reo%recv_count(sub_comm%size))
            allocate(reo%recv_offset(sub_comm%size))

            allocate(reo%dest(box%lo(dir):box%hi(dir)))


            !--------------------------------------------------------------
            ! Obtain 'reverted' bounds of all MPI ranks

            allocate(rlos(0:layout%size(dir)-1))
            allocate(rhis(0:layout%size(dir)-1))
            do i = 0, layout%size(dir)-1
                call get_local_bounds(ncell, i, layout%size(dir), lo, hi)
                rhis(i) = ncell - lo - 1
                rlos(i) = ncell - hi - 1
            enddo

            reo%send_recv_count = 0

            !--------------------------------------------------------------
            ! Determine destination rank and number of elements to send in x
            do i = box%lo(dir), box%hi(dir)
                reo%dest(i) = sub_comm%rank
                do j = 0, layout%size(dir)-1
                    if (i >= rlos(j) .and. i <= rhis(j)) then
!                         print *, sub_comm%rank, sub_comm%coord, j
                        call MPI_Cart_rank(sub_comm%comm, (/j/), rank, comm%err)
                        reo%dest(i) = rank
                    endif
                enddo
                j = reo%dest(i)
                if (.not. sub_comm%rank == j) then
                    reo%send_recv_count(j+1) = reo%send_recv_count(j+1) + 1
                endif
            enddo

            !--------------------------------------------------------------
            ! Add y and z dimensions
            reo%recv_count = reo%send_recv_count
            reo%send_recv_count = reo%send_recv_count * box%size(d) * box%size(I_Z)

            deallocate(rlos)
            deallocate(rhis)

            allocate(reo%send_buffer(sum(reo%send_recv_count)))
            allocate(reo%recv_buffer(sum(reo%send_recv_count)))

            !--------------------------------------------------------------
            ! Determine offsets for send and receive buffers

            reo%send_offset(1) = 0
            do i = 2, sub_comm%size
                reo%send_offset(i) = sum(reo%send_recv_count(1:i-1))
                reo%recv_offset(sub_comm%size-i+1) = reo%send_offset(i)
            enddo

            reo%recv_offset(sub_comm%size) = 0
            do i = sub_comm%size-1, 1, -1
                reo%recv_offset(i) = sum(reo%send_recv_count(sub_comm%size:i+1:-1))
            enddo

        end subroutine setup_reordering

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_x(fs)
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz, j

            j = 1
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%lo(2), box%hi(2)
                    do iz = box%lo(3), box%hi(3)
                        x_reo%send_buffer(j) = fs(iz, iy, ix)
                        j = j + 1
                    enddo
                enddo
            enddo

        end subroutine copy_to_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_x(fs)
            double precision, intent(out) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy, iz
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(3):box%hi(3), box%lo(2):box%hi(2))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(x_reo%recv_count)

            half_length = (box%size(1) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(1) + i
                k = box%hi(1) - n_recvs - i
                buf = fs(:, :, j)
                fs(:, :, j) = fs(:, :, k)
                fs(:, :, k) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(x_reo%recv_buffer)
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%hi(2), box%lo(2), -1
                    do iz = box%hi(3), box%lo(3), -1
                        fs(iz, iy, ix) = x_reo%recv_buffer(j)
                        j = j - 1
                    enddo
                enddo
            enddo

        end subroutine copy_from_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_y(fs)
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz, j

            j = 1
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(iy)) then
                    cycle
                endif
                do ix = box%lo(1), box%hi(1)
                    do iz = box%lo(3), box%hi(3)
                        y_reo%send_buffer(j) = fs(iz, iy, ix)
                        j = j + 1
                    enddo
                enddo
            enddo

        end subroutine copy_to_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_y(fs)
            double precision, intent(out) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy, iz
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(3):box%hi(3), box%lo(1):box%hi(1))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(y_reo%recv_count)

            half_length = (box%size(2) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(2) + i
                k = box%hi(2) - n_recvs - i
                buf = fs(:, j, :)
                fs(:, j, :) = fs(:, k, :)
                fs(:, k, :) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(y_reo%recv_buffer)
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(ix)) then
                    cycle
                endif
                do ix = box%hi(1), box%lo(1), -1
                    do iz = box%hi(3), box%lo(3), -1
                        fs(iz, iy, ix) = y_reo%recv_buffer(j)
                        j = j - 1
                    enddo
                enddo
            enddo

        end subroutine copy_from_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reorder(reo, sub_comm, fs, gs)
            type(reorder_type),                intent(inout) :: reo
            type(sub_communicator), intent(inout) :: sub_comm
            double precision,                  intent(in)    :: fs(box%hlo(3):box%hhi(3), &
                                                                   box%hlo(2):box%hhi(2), &
                                                                   box%hlo(1):box%hhi(1))
            double precision,                  intent(out)   :: gs(box%hlo(3):box%hhi(3), &
                                                                   box%hlo(2):box%hhi(2), &
                                                                   box%hlo(1):box%hhi(1))

            gs = fs

            call initialise_diff

            ! we exclude the halo when copying
            select case (reo%dir)
                case (I_X)
                    call copy_to_buffer_in_x(fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1)))
                case (I_Y)
                    call copy_to_buffer_in_y(fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1)))
                case default
                    ! FIXME
            end select

            call MPI_alltoallv(reo%send_buffer,         &
                               reo%send_recv_count,     &
                               reo%send_offset,         &
                               MPI_DOUBLE_PRECISION,    &
                               reo%recv_buffer,         &
                               reo%send_recv_count,     &
                               reo%recv_offset,         &
                               MPI_DOUBLE_PRECISION,    &
                               sub_comm%comm,           &
                               sub_comm%err)

            ! we exclude the halo when copying
            select case (reo%dir)
                case (I_X)
                    call copy_from_buffer_in_x(gs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1)))
                case (I_Y)
                    call copy_from_buffer_in_y(gs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1)))
                case default
                    ! FIXME
            end select

            call field_halo_fill(gs)

        end subroutine reorder

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs, ds)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: ds(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision              :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            integer                       :: kx, dkx, kxc, nwx, nxp2
            double precision              :: si

            nwx = nx / 2
            nxp2 = nx + 1

            call reorder(x_reo, x_comm, fs, gs)

            !Carry out differentiation by wavenumber multiplication:
            if (0 == box%lo(1)) then
                ds(:, :, 0) = zero
            endif

            do kx = max(1, box%lo(1)), box%hi(1)
                dkx = min(2 * kx, 2 * (nx - kx))
                si = merge(1.0d0, -1.0d0, kx >= nwx + 1)
                ds(:, :, kx)  = si * hrkx(dkx) * gs(:, :, kx-1)
            enddo

            if (mod(nx, 2) .eq. 0) then
                kxc = nwx! + 1
                if (kxc >= box%lo(1) .and. kxc <= box%hi(1)) then
                    ds(:, :, kxc) = zero
                endif
            endif

        end subroutine

end module fft_utils
