module fft_utils
    use mpi_communicator
    use mpi_layout
    use constants, only : zero
    use field_mpi, only : field_halo_fill
    use parameters, only : nx
    use inversion_utils, only : hrkx
    implicit none

    integer, allocatable :: dest(:)                           ! destination rank
    integer, allocatable :: send_count(:), send_offset(:)
    integer, allocatable :: recv_count(:), recv_offset(:), x_recv_count(:)
    double precision, allocatable :: send_buffer(:), recv_buffer(:)

    logical :: l_initialised = .false.

    private :: copy_from_buffer, copy_to_buffer, l_initialised

    type(MPI_type) :: comm_test
    integer :: rr

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine initialise_fft_diff
            integer, allocatable :: rlos(:), rhis(:)
            integer              :: i, j, lo, hi, rank

            if (l_initialised) then
                return
            endif
            l_initialised = .true.

            allocate(send_count(comm%size))
            allocate(send_offset(comm%size))
            allocate(recv_count(comm%size))
            allocate(x_recv_count(comm%size))
            allocate(recv_offset(comm%size))

            allocate(dest(box%lo(1):box%hi(1)))


            !--------------------------------------------------------------
            ! Obtain 'reverted' bounds of all MPI ranks

            allocate(rlos(0:layout%size(1)-1))
            allocate(rhis(0:layout%size(1)-1))
            do i = 0, layout%size(1)-1
                call get_local_bounds(nx, i, layout%size(1), lo, hi)
                rhis(i) = nx - lo - 1
                rlos(i) = nx - hi - 1
!                 print *, comm%rank, "lo:", lo, "hi:", hi
            enddo

!             print *, comm%rank, "my:", box%lo(1), box%hi(1)
!             print *, comm%rank, "re:", rhis(comm%rank), rlos(comm%rank)

            send_count(:) = 0

            !--------------------------------------------------------------
            ! Determine destination rank and number of elements to send in x
            do i = box%lo(1), box%hi(1)
                dest(i) = comm%rank
                do j = 0, layout%size(1)-1
                    if (i >= rlos(j) .and. i <= rhis(j)) then
                        call MPI_Cart_rank(comm%cart, (/j, layout%coords(2)/), rank, comm%err)
                        dest(i) = rank
                    endif
                enddo
                j = dest(i)
                if (.not. comm%rank == j) then
                    send_count(j+1) = send_count(j+1) + 1
                endif
            enddo

            cat MPI_Cart_sub(comm%cart, (/.true., .false./), comm_test, comm%err)

            call MPI_Comm_rank(comm_test, rr, comm%err)
            print *, comm%rank, rr
!             call MPI_Cart_coords(comm_test, rank, maxdims, coords, ierror)

            print *, comm%rank, "dest", dest

            !--------------------------------------------------------------
            ! Add y and z dimensions
            x_recv_count = send_count
            send_count = send_count * box%size(2) * box%size(3)
            recv_count = send_count

            deallocate(rlos)
            deallocate(rhis)

            allocate(send_buffer(sum(send_count)))
            allocate(recv_buffer(sum(recv_count)))

            !--------------------------------------------------------------
            ! Determine offsets for send and receive buffers

            send_offset(1) = 0
            do i = 2, comm%size
                send_offset(i) = sum(send_count(1:i-1))
                recv_offset(comm%size-i+1) = send_offset(i)
            enddo

            recv_offset(comm%size) = 0
            do i = comm%size-1, 1, -1
                recv_offset(i) = sum(recv_count(comm%size:i+1:-1))
            enddo

            print *, comm%rank, "send_count", send_count
            print *, comm%rank, "recv_count", recv_count
            print *, comm%rank, "send_offset", send_offset
            print *, comm%rank, "recv_offset", recv_offset

        end subroutine initialise_fft_diff

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer(fs)
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz, j

            j = 1
            do ix = box%hi(1), box%lo(1), -1
                if (comm%rank == dest(ix)) then
                    cycle
                endif
                do iy = box%lo(2), box%hi(2)
                    do iz = box%lo(3), box%hi(3)
                        send_buffer(j) = fs(iz, iy, ix)
                        j = j + 1
                    enddo
                enddo
            enddo

        end subroutine copy_to_buffer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer(fs)
            double precision, intent(out) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy, iz
            integer                       :: i, j, k, half_length, n_recvs
            double precision, allocatable :: buf(:, :)


            !--------------------------------------------------------------
            ! Revert array locally

            allocate(buf(box%lo(3):box%hi(3), box%lo(2):box%hi(2)))

            n_recvs = sum(x_recv_count)

            print *, comm%rank, "n_recvs", n_recvs

            half_length = (box%size(1) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(1) + i
                k = box%hi(1) - n_recvs - i
                buf = fs(:, :, j)
                fs(:, :, j) = fs(:, :, k)
                fs(:, :, k) = buf
            enddo

            deallocate(buf)

            print *, comm%rank, "hs:", fs(0, 0, :)

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(recv_buffer)
            do ix = box%hi(1), box%lo(1), -1
                if (comm%rank == dest(ix)) then
                    cycle
                endif
                do iy = box%hi(2), box%lo(2), -1
                    do iz = box%hi(3), box%lo(3), -1
                        fs(iz, iy, ix) = recv_buffer(j)
                        j = j - 1
                    enddo
                enddo
            enddo

            print *, comm%rank, "fs:", fs(0, 0, :)

        end subroutine copy_from_buffer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reorder(fs, gs)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))

            gs = fs

            call initialise_fft_diff

            ! we exclude the halo when copying
            call copy_to_buffer(fs(box%lo(3):box%hi(3), &
                                   box%lo(2):box%hi(2), &
                                   box%lo(1):box%hi(1)))

            call MPI_alltoallv(send_buffer,             &
                               send_count,              &
                               send_offset,             &
                               MPI_DOUBLE_PRECISION,    &
                               recv_buffer,             &
                               recv_count,              &
                               recv_offset,             &
                               MPI_DOUBLE_PRECISION,    &
                               comm%cart,               &
                               comm%err)

            ! we exclude the halo when copying
            call copy_from_buffer(gs(box%lo(3):box%hi(3), &
                                     box%lo(2):box%hi(2), &
                                     box%lo(1):box%hi(1)))

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

            call reorder(fs, gs)

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
