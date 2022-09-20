module parcel_mpi
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error
    use mpi_tags
    use fields, only : get_index
    use merge_sort, only : msort
    use parameters, only : vmin, vcell, nz
    use parcel_container, only : n_par_attrib       &
                               , parcel_serialize   &
                               , parcel_deserialize &
                               , parcel_delete      &
                               , n_parcels          &
                               , parcels
    implicit none

    private

    integer, allocatable, dimension(:) :: north_buf, south_buf, west_buf, east_buf, &
                                          northwest_buf, northeast_buf,             &
                                          southwest_buf, southeast_buf,             &
                                          invalid

    ! we have 8 neighbours
    integer :: n_sends(8)
    integer :: n_recvs(8)

    public :: parcel_halo_swap  &
            , pack_parcels      &
            , unpack_parcels    &
            , north_buf         &
            , south_buf         &
            , west_buf          &
            , east_buf          &
            , northwest_buf     &
            , northeast_buf     &
            , southwest_buf     &
            , southeast_buf

    contains

        subroutine parcel_halo_swap(pindex)
            integer, optional, intent(in) :: pindex(:)
            integer                       :: n_total_sends
            type(MPI_Request)             :: request

            call allocate_buffers

            ! figure out where parcels should go
            call locate_parcels(pindex)

            ! tell your neighbours the number of receiving parcels
            call MPI_Isend(n_sends(NB_NORTH), 1, MPI_INT, neighbour%north, NORTH_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_SOUTH), 1, MPI_INT, neighbour%south, NORTH_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_SOUTH), 1, MPI_INT, neighbour%south, SOUTH_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_NORTH), 1, MPI_INT, neighbour%north, SOUTH_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_WEST), 1, MPI_INT, neighbour%west, WEST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_EAST), 1, MPI_INT, neighbour%east, WEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_EAST), 1, MPI_INT, neighbour%east, EAST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_WEST), 1, MPI_INT, neighbour%west, EAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, NORTHWEST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, NORTHWEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, SOUTHEAST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, SOUTHEAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, NORTHEAST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, NORTHEAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Isend(n_sends(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, SOUTHWEST_TAG, &
                           comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(n_recvs(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, SOUTHWEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            ! communicate parcels
            call exchange_parcels(north_buf, n_sends(NB_NORTH), neighbour%north, &
                                             n_recvs(NB_SOUTH), neighbour%south, NORTH_TAG)

            call exchange_parcels(south_buf, n_sends(NB_SOUTH), neighbour%south, &
                                             n_recvs(NB_NORTH), neighbour%north, SOUTH_TAG)

            call exchange_parcels(west_buf, n_sends(NB_WEST), neighbour%west, &
                                            n_recvs(NB_EAST), neighbour%east, WEST_TAG)

            call exchange_parcels(east_buf, n_sends(NB_EAST), neighbour%east, &
                                            n_recvs(NB_WEST), neighbour%west, EAST_TAG)

            call exchange_parcels(northwest_buf, &
                                  n_sends(NB_NORTHWEST), neighbour%northwest, &
                                  n_recvs(NB_SOUTHEAST), neighbour%southeast, NORTHWEST_TAG)

            call exchange_parcels(northeast_buf, &
                                  n_sends(NB_NORTHEAST), neighbour%northeast, &
                                  n_recvs(NB_SOUTHWEST), neighbour%southwest, SOUTHEAST_TAG)

            call exchange_parcels(southwest_buf,  &
                                  n_sends(NB_SOUTHWEST), neighbour%southwest, &
                                  n_recvs(NB_NORTHEAST), neighbour%northeast, NORTHEAST_TAG)

            call exchange_parcels(southeast_buf, &
                                  n_sends(NB_SOUTHEAST), neighbour%southeast, &
                                  n_recvs(NB_NORTHWEST), neighbour%northwest, SOUTHWEST_TAG)

            ! delete parcel that we sent
            n_total_sends = sum(n_sends)
            call parcel_delete(invalid, n_total_sends)

        end subroutine parcel_halo_swap


        subroutine exchange_parcels(pid, sendcount, dest, recvcount, source, tag)
            integer, intent(in)           :: pid(:)
            integer, intent(in)           :: sendcount, recvcount
            integer, intent(in)           :: tag
            integer, intent(in)           :: dest, source
            double precision, allocatable :: sendbuf(:), recvbuf(:)
            integer                       :: send_size, recv_size
            type(MPI_Request)             :: request

            send_size = n_par_attrib * sendcount
            recv_size = n_par_attrib * recvcount

            allocate(sendbuf(send_size))
            allocate(recvbuf(recv_size))

            call pack_parcels(pid, sendcount, sendbuf)

            call MPI_Isend(sendbuf, send_size, MPI_DOUBLE, dest, tag, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(recvbuf, recv_size, MPI_DOUBLE, source, tag, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call unpack_parcels(recvcount, recvbuf)

            deallocate(sendbuf)
            deallocate(recvbuf)

        end subroutine exchange_parcels


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
                    case (NB_NONE)
                        ! do nothing
                        cycle
                    case (NB_NORTH)
                        m = n_sends(nb) + 1
                        north_buf(m) = l
                    case (NB_SOUTH)
                        m = n_sends(nb) + 1
                        south_buf(m) = l
                    case (NB_WEST)
                        m = n_sends(nb) + 1
                        west_buf(m) = l
                    case (NB_EAST)
                        m = n_sends(nb) + 1
                        east_buf(m) = l
                    case (NB_NORTHWEST)
                        m = n_sends(nb) + 1
                        northwest_buf(m) = l
                    case (NB_NORTHEAST)
                        m = n_sends(nb) + 1
                        northeast_buf(m) = l
                    case (NB_SOUTHWEST)
                        m = n_sends(nb) + 1
                        southwest_buf(m) = l
                    case (NB_SOUTHEAST)
                        m = n_sends(nb) + 1
                        southeast_buf(m) = l
                    case default
                        call mpi_exit_on_error("locate_parcels: Parcel was not assigned properly.")
                end select

                invalid(iv) = l
                n_sends(nb) = n_sends(nb) + 1
                iv = iv + 1
            enddo

        end subroutine locate_parcels


        subroutine pack_parcels(pid, sendcount, sendbuf)
            integer,          intent(in)  :: pid(:)
            integer,          intent(in)  :: sendcount
            double precision, intent(out) :: sendbuf(:)
            integer                       :: n, i, j

            do n = 1, sendcount
                i = 1 + (n-1) * n_par_attrib
                j = n * n_par_attrib
                call parcel_serialize(pid(n), sendbuf(i:j))
            enddo
        end subroutine pack_parcels


        subroutine unpack_parcels(recvcount, recvbuf)
            integer,          intent(in) :: recvcount
            double precision, intent(in) :: recvbuf(:)
            integer                      :: n, i, j

            do n = 1, recvcount
                i = 1 + (n-1) * n_par_attrib
                j = n * n_par_attrib
                call parcel_deserialize(n_parcels + n, recvbuf(i:j))
            enddo

            n_parcels = n_parcels + recvcount

        end subroutine unpack_parcels


        subroutine allocate_buffers
            integer :: nc, ub, n_max

            if (allocated(north_buf)) then
                return
            endif

            ! upper bound (ub) of parcels per cell
            ub = ceiling(vcell / vmin)

            ! number of cells sharing with north and south neighbour
            nc = (box%hi(1) - box%lo(1) + 1) * nz
            n_max = 2 * nc

            allocate(north_buf(nc * ub))
            allocate(south_buf(nc * ub))

            ! number of cells sharing with west and east neighbour
            nc = (box%hi(2) - box%lo(2) + 1) * nz
            n_max = n_max + 2 * nc

            allocate(west_buf(nc * ub))
            allocate(east_buf(nc * ub))

            ! number of cells sharing with corner neighbours
            nc = nz
            n_max = n_max + 4 * nc

            allocate(northwest_buf(nc * ub))
            allocate(northeast_buf(nc * ub))
            allocate(southwest_buf(nc * ub))
            allocate(southeast_buf(nc * ub))

            ! Note: This buffer needs to have start index 0 because of the parcel delete routine.
            !       However, this first array element is never assigned any value.
            allocate(invalid(0:n_max * ub))

        end subroutine allocate_buffers

end module parcel_mpi
