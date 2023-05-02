module parcel_mpi
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error     &
                        , mpi_check_for_message &
                        , mpi_check_for_error
    use mpi_tags
#ifndef NDEBUG
    use mpi_collectives, only : mpi_blocking_reduce
#endif
    use fields, only : get_index_periodic
    use merge_sort, only : msort
    use parameters, only : vmin, vcell, nz
    use parcel_container, only : n_par_attrib       &
                               , parcel_pack        &
                               , parcel_unpack      &
                               , parcel_delete      &
                               , n_parcels          &
#ifndef NDEBUG
                               , n_total_parcels    &
#endif
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
    integer :: n_parcel_sends(8)

    public :: parcel_communicate,             &
              n_parcel_sends,               &
              north_pid,                    &
              south_pid,                    &
              west_pid,                     &
              east_pid,                     &
              northwest_pid,                &
              northeast_pid,                &
              southwest_pid,                &
              southeast_pid,                &
              north_buf,                    &
              south_buf,                    &
              west_buf,                     &
              east_buf,                     &
              northwest_buf,                &
              northeast_buf,                &
              southwest_buf,                &
              southeast_buf,                &
              invalid,                      &
              allocate_parcel_buffers,      &
              deallocate_parcel_buffers,    &
              allocate_parcel_id_buffers,   &
              deallocate_parcel_id_buffers, &
              get_parcel_buffer_ptr,        &
              communicate_parcels

    contains

        ! We follow here PMPIC
        ! https://github.com/EPCCed/pmpic/blob/master/model_core/src/parcels/haloswap.F90
        ! git sha: 4a77642ef82096a770324f5f1dc587ce121065ab
        subroutine get_parcel_buffer_ptr(dir, pid_ptr, buf_ptr)
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
                    call mpi_exit_on_error(&
                        "in parcel_mpi::get_parcel_buffer_ptr: No valid direction.")
            end select

        end subroutine get_parcel_buffer_ptr

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_communicate(pindex)
            integer, optional, intent(in)           :: pindex(:)

            ! We only must store the parcel indices (therefore 1) and
            ! also allocate the buffer for invalid parcels. (therefore .true.)
            call allocate_parcel_id_buffers(1, .true.)

            ! figure out where parcels go
            call locate_parcels(pindex)

            call allocate_parcel_buffers(n_par_attrib)

            call communicate_parcels

            call deallocate_parcel_id_buffers

            call deallocate_parcel_buffers

        end subroutine parcel_communicate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine communicate_parcels
            integer,          dimension(:), pointer :: pid
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: n_total_sends, n
            integer                                 :: recv_count
            integer                                 :: recv_size, send_size

            do n = 1, 8
                call get_parcel_buffer_ptr(n, pid, send_buf)

                send_size = n_parcel_sends(n) * n_par_attrib

                if (n_parcel_sends(n) > 0) then
                    call parcel_pack(pid, n_parcel_sends(n), send_buf)
                endif

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error("in MPI_Isend of parcel_mpi::communicate_parcels.")
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                call mpi_check_for_error("in MPI_Recv of parcel_mpi::communicate_parcels.")

                if (mod(recv_size, n_par_attrib) /= 0) then
                    call mpi_exit_on_error("parcel_mpi::communicate_parcels: Receiving wrong count.")
                endif

                recv_count = recv_size / n_par_attrib

                if (recv_count > 0) then
                    call parcel_unpack(recv_count, recv_buf)
                endif

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            call mpi_check_for_error("in MPI_Waitall of parcel_mpi::communicate_parcels.")

            ! delete parcel that we sent
            n_total_sends = sum(n_parcel_sends)
            call parcel_delete(invalid, n_total_sends)

#ifndef NDEBUG
            n = n_parcels
            call mpi_blocking_reduce(n, MPI_SUM)
            if ((comm%rank == comm%master) .and. (.not. n == n_total_parcels)) then
                call mpi_exit_on_error(&
                    "in parcel_mpi::communicate_parcels: We lost parcels.")
            endif
#endif
        end subroutine communicate_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine locate_parcels(pindex)
            integer, optional, intent(in) :: pindex(:)
            integer                       :: i, j, k, n, nb, m, iv, np, l

            ! reset the number of sends
            n_parcel_sends(:) = 0

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

                call get_index_periodic(parcels%position(:, l), i, j, k)

                nb = get_neighbour(i, j)

                select case (nb)
                    case (MPI_NONE)
                        ! do nothing
                        cycle
                    case (MPI_NORTH)
                        m = n_parcel_sends(nb) + 1
                        north_pid(m) = l
                    case (MPI_SOUTH)
                        m = n_parcel_sends(nb) + 1
                        south_pid(m) = l
                    case (MPI_WEST)
                        m = n_parcel_sends(nb) + 1
                        west_pid(m) = l
                    case (MPI_EAST)
                        m = n_parcel_sends(nb) + 1
                        east_pid(m) = l
                    case (MPI_NORTHWEST)
                        m = n_parcel_sends(nb) + 1
                        northwest_pid(m) = l
                    case (MPI_NORTHEAST)
                        m = n_parcel_sends(nb) + 1
                        northeast_pid(m) = l
                    case (MPI_SOUTHWEST)
                        m = n_parcel_sends(nb) + 1
                        southwest_pid(m) = l
                    case (MPI_SOUTHEAST)
                        m = n_parcel_sends(nb) + 1
                        southeast_pid(m) = l
                    case default
                        call mpi_exit_on_error("locate_parcels: Parcel was not assigned properly.")
                end select

                invalid(iv) = l
                n_parcel_sends(nb) = n_parcel_sends(nb) + 1
                iv = iv + 1
            enddo

        end subroutine locate_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_parcel_id_buffers(n_ids, l_invalid)
            integer, intent(in) :: n_ids
            logical, intent(in) :: l_invalid
            integer             :: nc, ub, n_max

            n_parcel_sends = 0

            ! upper bound (ub) of parcels per cell
            ub = ceiling(vcell / vmin)

            ! number of cells sharing with north and south neighbour
            nc = (box%hi(1) - box%lo(1) + 1) * nz
            n_max = 2 * nc

            allocate(north_pid(nc * ub * n_ids))
            allocate(south_pid(nc * ub * n_ids))

            ! number of cells sharing with west and east neighbour
            nc = (box%hi(2) - box%lo(2) + 1) * nz
            n_max = n_max + 2 * nc

            allocate(west_pid(nc * ub * n_ids))
            allocate(east_pid(nc * ub * n_ids))

            ! number of cells sharing with corner neighbours
            nc = nz
            n_max = n_max + 4 * nc

            allocate(northwest_pid(nc * ub * n_ids))
            allocate(northeast_pid(nc * ub * n_ids))
            allocate(southwest_pid(nc * ub * n_ids))
            allocate(southeast_pid(nc * ub * n_ids))

            if (l_invalid) then
                ! Note: This buffer needs to have start index 0 because of the parcel delete routine.
                !       However, this first array element is never assigned any value.
                allocate(invalid(0:n_max * ub * n_ids))
            endif

        end subroutine allocate_parcel_id_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine deallocate_parcel_id_buffers
            deallocate(north_pid)
            deallocate(south_pid)
            deallocate(west_pid)
            deallocate(east_pid)

            deallocate(northwest_pid)
            deallocate(northeast_pid)
            deallocate(southwest_pid)
            deallocate(southeast_pid)

            if (allocated(invalid)) then
                deallocate(invalid)
            endif
        end subroutine deallocate_parcel_id_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_parcel_buffers(n_attributes)
            integer, intent(in) :: n_attributes

            allocate(north_buf(n_parcel_sends(MPI_NORTH) * n_attributes))
            allocate(northeast_buf(n_parcel_sends(MPI_NORTHEAST) * n_attributes))
            allocate(east_buf(n_parcel_sends(MPI_EAST) * n_attributes))
            allocate(southeast_buf(n_parcel_sends(MPI_SOUTHEAST) * n_attributes))
            allocate(south_buf(n_parcel_sends(MPI_SOUTH) * n_attributes))
            allocate(southwest_buf(n_parcel_sends(MPI_SOUTHWEST) * n_attributes))
            allocate(west_buf(n_parcel_sends(MPI_WEST) * n_attributes))
            allocate(northwest_buf(n_parcel_sends(MPI_NORTHWEST) * n_attributes))
        end subroutine allocate_parcel_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine deallocate_parcel_buffers
            deallocate(north_buf)
            deallocate(northeast_buf)
            deallocate(east_buf)
            deallocate(southeast_buf)
            deallocate(south_buf)
            deallocate(southwest_buf)
            deallocate(west_buf)
            deallocate(northwest_buf)
        end subroutine deallocate_parcel_buffers

end module parcel_mpi
