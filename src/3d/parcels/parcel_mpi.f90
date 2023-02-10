module parcel_mpi
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_message
    use mpi_tags
    use fields, only : get_index
    use merge_sort, only : msort
    use parameters, only : vmin, vcell, nz
    use parcel_container, only : n_par_attrib       &
                               , parcel_pack        &
                               , parcel_unpack      &
                               , parcel_delete      &
                               , n_parcels          &
                               , parcels
    implicit none

    private

    integer, allocatable, dimension(:), target :: north_pid, south_pid, west_pid, east_pid, &
                                                  northwest_pid, northeast_pid,             &
                                                  southwest_pid, southeast_pid
    double precision, allocatable, dimension(:), target :: north_buf, south_buf, west_buf, east_buf, &
                                                           northwest_buf, northeast_buf,             &
                                                           southwest_buf, southeast_buf
    integer, allocatable, dimension(:) :: invalid

    ! we have 8 neighbours
    integer :: n_sends(8)

    public :: parcel_halo_swap

    contains

        ! We follow here PMPIC
        ! https://github.com/EPCCed/pmpic/blob/master/model_core/src/parcels/haloswap.F90
        ! git sha: 4a77642ef82096a770324f5f1dc587ce121065ab
        subroutine get_buffer_ptr(dir, pid_ptr, buf_ptr)
            integer,                                 intent(in)  :: dir
            integer,          dimension(:), pointer, intent(out) :: pid_ptr
            double precision, dimension(:), pointer, intent(out) :: buf_ptr

            select case (dir)
                case (MPI_NORTH)
                    pid_ptr => north_pid
                    buf_ptr => north_buf
                case (MPI_SOUTH)
                    pid_ptr => south_pid
                    buf_ptr => south_buf
                case (MPI_WEST)
                    pid_ptr => west_pid
                    buf_ptr => west_buf
                case (MPI_EAST)
                    pid_ptr => east_pid
                    buf_ptr => east_buf
                case (MPI_NORTHWEST)
                    pid_ptr => northwest_pid
                    buf_ptr => northwest_buf
                case (MPI_NORTHEAST)
                    pid_ptr => northeast_pid
                    buf_ptr => northeast_buf
                case (MPI_SOUTHWEST)
                    pid_ptr => southwest_pid
                    buf_ptr => southwest_buf
                case (MPI_SOUTHEAST)
                    pid_ptr => southeast_pid
                    buf_ptr => southeast_buf
                case default
                    call mpi_exit_on_error("get_buffer_ptr: No valid direction.")
            end select

        end subroutine get_buffer_ptr

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_halo_swap(pindex)
            integer, optional, intent(in)           :: pindex(:)
            integer,          dimension(:), pointer :: pid
            double precision, dimension(:), pointer :: sendbuf
            double precision, allocatable           :: recvbuf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: n_total_sends, n
            integer                                 :: tag, source, recvcount
            integer                                 :: recv_size, send_size

            call allocate_pid_buffers

            ! figure out where parcels go
            call locate_parcels(pindex)

            call allocate_send_buffers

            do n = 1, 8
                call get_buffer_ptr(n, pid, sendbuf)

                send_size = n_sends(n) * n_par_attrib

                if (n_sends(n) > 0) then
                    call parcel_pack(pid, n_sends(n), sendbuf)
                endif

                call MPI_Isend(sendbuf,                 &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                allocate(recvbuf(recv_size))

                call MPI_Recv(recvbuf,                  &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                recvcount = recv_size / n_par_attrib

                if (recvcount > 0) then
                    call parcel_unpack(recvcount, recvbuf)
                endif

                deallocate(recvbuf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            ! delete parcel that we sent
            n_total_sends = sum(n_sends)
            call parcel_delete(invalid, n_total_sends)

            call deallocate_pid_buffers

            call deallocate_send_buffers

        end subroutine parcel_halo_swap

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine locate_parcels(pindex)
            integer, optional, intent(in) :: pindex(:)
            integer                       :: i, j, k, n, nb, m, iv, np, l

            ! reset the number of sends
            n_sends(:) = 0

            np = n_parcels

            if (present(pindex)) then
                np = size(pindex)
            endif

            iv = 1
            do n = 1, np

                l = n
                if (present(pindex)) then
                    l = pindex(n)
                endif

                call get_index(parcels%position(:, l), i, j, k)

                nb = get_neighbour(i, j)

                select case (nb)
                    case (MPI_NONE)
                        ! do nothing
                        cycle
                    case (MPI_NORTH)
                        m = n_sends(nb) + 1
                        north_pid(m) = l
                    case (MPI_SOUTH)
                        m = n_sends(nb) + 1
                        south_pid(m) = l
                    case (MPI_WEST)
                        m = n_sends(nb) + 1
                        west_pid(m) = l
                    case (MPI_EAST)
                        m = n_sends(nb) + 1
                        east_pid(m) = l
                    case (MPI_NORTHWEST)
                        m = n_sends(nb) + 1
                        northwest_pid(m) = l
                    case (MPI_NORTHEAST)
                        m = n_sends(nb) + 1
                        northeast_pid(m) = l
                    case (MPI_SOUTHWEST)
                        m = n_sends(nb) + 1
                        southwest_pid(m) = l
                    case (MPI_SOUTHEAST)
                        m = n_sends(nb) + 1
                        southeast_pid(m) = l
                    case default
                        call mpi_exit_on_error("locate_parcels: Parcel was not assigned properly.")
                end select

                invalid(iv) = l
                n_sends(nb) = n_sends(nb) + 1
                iv = iv + 1
            enddo

        end subroutine locate_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_pid_buffers
            integer :: nc, ub, n_max

            ! upper bound (ub) of parcels per cell
            ub = ceiling(vcell / vmin)

            ! number of cells sharing with north and south neighbour
            nc = (box%hi(1) - box%lo(1) + 1) * nz
            n_max = 2 * nc

            allocate(north_pid(nc * ub))
            allocate(south_pid(nc * ub))

            ! number of cells sharing with west and east neighbour
            nc = (box%hi(2) - box%lo(2) + 1) * nz
            n_max = n_max + 2 * nc

            allocate(west_pid(nc * ub))
            allocate(east_pid(nc * ub))

            ! number of cells sharing with corner neighbours
            nc = nz
            n_max = n_max + 4 * nc

            allocate(northwest_pid(nc * ub))
            allocate(northeast_pid(nc * ub))
            allocate(southwest_pid(nc * ub))
            allocate(southeast_pid(nc * ub))

            ! Note: This buffer needs to have start index 0 because of the parcel delete routine.
            !       However, this first array element is never assigned any value.
            allocate(invalid(0:n_max * ub))

        end subroutine allocate_pid_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine deallocate_pid_buffers
            deallocate(north_pid)
            deallocate(south_pid)
            deallocate(west_pid)
            deallocate(east_pid)

            deallocate(northwest_pid)
            deallocate(northeast_pid)
            deallocate(southwest_pid)
            deallocate(southeast_pid)

            deallocate(invalid)
        end subroutine deallocate_pid_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_send_buffers
            allocate(north_buf(n_sends(MPI_NORTH) * n_par_attrib))
            allocate(northeast_buf(n_sends(MPI_NORTHEAST) * n_par_attrib))
            allocate(east_buf(n_sends(MPI_EAST) * n_par_attrib))
            allocate(southeast_buf(n_sends(MPI_SOUTHEAST) * n_par_attrib))
            allocate(south_buf(n_sends(MPI_SOUTH) * n_par_attrib))
            allocate(southwest_buf(n_sends(MPI_SOUTHWEST) * n_par_attrib))
            allocate(west_buf(n_sends(MPI_WEST) * n_par_attrib))
            allocate(northwest_buf(n_sends(MPI_NORTHWEST) * n_par_attrib))
        end subroutine allocate_send_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine deallocate_send_buffers
            deallocate(north_buf)
            deallocate(northeast_buf)
            deallocate(east_buf)
            deallocate(southeast_buf)
            deallocate(south_buf)
            deallocate(southwest_buf)
            deallocate(west_buf)
            deallocate(northwest_buf)
        end subroutine deallocate_send_buffers

end module parcel_mpi
