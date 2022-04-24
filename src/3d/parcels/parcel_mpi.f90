module parcel_mpi
    use mpi_layout, only : box
    use mpi_utils, only : mpi_exit_on_error
    use mpi_tags
    use parcel_container, only : n_par_attrib, parcel_serialize, parcel_deserialize, n_parcels
    implicit none

    ! buffers storing the parcel indices the neighbours get
    integer, allocatable :: northwest(:), north(:), northeast(:), west(:) &
                            southwest(:), south(:), southeast(:), east(:)

    ! we have 8 neighbours (the index 0, means the parcels is kept)
    integer :: n_sends(0:8)
    integer :: n_recvs(8)

    ! current last index in parcel container
    integer :: n_last

    public :: parcel_halo_swap

    contains

        subroutine parcel_halo_swap

            call allocate_buffers

            ! figure out where parcels should go
            call locate_parcels

            ! tell your neighbours the number of receiving parcels
            call MPI_Sendrecv(n_sends(NB_NORTH), 1, MPI_INT, neighbour%north, SEND_NORTH_TAG, &
                              n_recvs(NB_SOUTH), 1, MPI_INT, neighbour%south, RECV_SOUTH_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_SOUTH), 1, MPI_INT, neighbour%south, SEND_SOUTH_TAG, &
                              n_recvs(NB_NORTH), 1, MPI_INT, neighbour%north, RECV_NORTH_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_WEST), 1, MPI_INT, neighbour%west, SEND_WEST_TAG, &
                              n_recvs(NB_EAST), 1, MPI_INT, neighbour%east, RECV_EAST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_EAST), 1, MPI_INT, neighbour%east, SEND_EAST_TAG, &
                              n_recvs(NB_WEST), 1, MPI_INT, neighbour%west, RECV_WEST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, SEND_NORTHWEST_TAG, &
                              n_recvs(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, RECV_SOUTHEAST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, SEND_SOUTHEAST_TAG, &
                              n_recvs(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, RECV_NORTHWEST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, SEND_NORTHEAST_TAG, &
                              n_recvs(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, RECV_SOUTHWEST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Sendrecv(n_sends(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, SEND_SOUTHWEST_TAG, &
                              n_recvs(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, RECV_NORTHEAST_TAG, &
                              comm_cart, MPI_STATUS_IGNORE, mpi_err)

            n_last = n_parcels

            ! communicate parcels
            call exchange_parcels(north, n_sends(NB_NORTH), neighbour%north, SEND_NORTH_TAG, &
                                         n_recvs(NB_SOUTH), neighbour%south, RECV_SOUTH_TAG)

            call exchange_parcels(south, n_sends(NB_SOUTH), neighbour%south, SEND_SOUTH_TAG, &
                                         n_recvs(NB_NORTH), neighbour%north, RECV_NORTH_TAG)

            call exchange_parcels(west, n_sends(NB_WEST), neighbour%west, SEND_WEST_TAG, &
                                        n_recvs(NB_EAST), neighbour%east, RECV_EAST_TAG)

            call exchange_parcels(east, n_sends(NB_EAST), neighbour%east, SEND_EAST_TAG, &
                                        n_recvs(NB_WEST), neighbour%west, RECV_WEST_TAG)

            call exchange_parcels(northwest, n_sends(NB_NORTHWEST), neighbour%northwest, SEND_NORTHWEST_TAG, &
                                             n_recvs(NB_SOUTHEAST), neighbour%southeast, RECV_SOUTHEAST_TAG)

            call exchange_parcels(northeast, n_sends(NB_NORTHEAST), neighbour%northeast, SEND_SOUTHEAST_TAG, &
                                             n_recvs(NB_SOUTHWEST), neighbour%southwest, RECV_NORTHWEST_TAG)

            call exchange_parcels(southwest, n_sends(NB_SOUTHWEST), neighbour%southwest, SEND_NORTHEAST_TAG, &
                                             n_recvs(NB_NORTHEAST), neighbour%northeast, RECV_SOUTHWEST_TAG)

            call exchange_parcels(southeast, n_sends(NB_SOUTHEAST), neighbour%southeast, SEND_SOUTHWEST_TAG, &
                                             n_recvs(NB_NORTHWEST), neighbour%northwest, RECV_NORTHEAST_TAG)


            ! delete parcels that we sent

        end subroutine parcel_halo_swap


        subroutine exchange_parcels(pid, sendcount, dest, sendtag, recvcount, source, recvtag)
            integer, intent(in)           :: pid(:)
            integer, intent(in)           :: sendcount, recvcount
            integer, intent(in)           :: sendtag, recvtag
            integer, intent(in)           :: dest, source
            double precision, allocatable :: sendbuf(:), recvbuf(:)
            integer                       :: send_size, recv_size

            send_size = n_par_attrib * sendcount
            recv_size = n_par_attrib * recvcount

            allocate(sendbuf(send_size))
            allocate(recvbuf(recv_size))

            call pack_parcels(pid, sendcount, sendbuf)

            call MPI_Sendrecv(sendbuf   = sendbuf,              &
                              sendcount = send_size,            &
                              sendtype  = MPI_DOUBLE,           &
                              dest      = dest,                 &
                              sendtag   = sendtag,              &
                              recvbuf   = recbuf,               &
                              recvcount = recv_size,            &
                              recvtype  = MPI_DOUBLE,           &
                              source    = source,               &
                              recvtag   = recvtag,              &
                              comm      = comm_cart,            &
                              status    = MPI_IGNORE_STATUS,    &
                              ierror    = mpi_err)

            call unpack_parcels(recvcount, recvbuf)

            deallocate(sendbuf)
            deallocate(recvbuf)

        end subroutine exchange_parcels


        subroutine locate_parcels
            integer              :: i, j, k, n, nb, idx

            ! reset the number of sends
            n_sends = 0

            do n = 1, n_parcels
                call get_index(parcels%position(:, n), i, j, k)

                nb = get_neighbour(i, j)

                idx = n_sends(nb)

                select case (nb)
                    case (NB_NONE)
                        ! parcel still owned by this rank
                    case (NB_NORTH)
                        north(idx) = n
                    case (NB_SOUTH)
                        south(idx) = n
                    case (NB_WEST)
                        west(idx) = n
                    case (NB_EAST)
                        east(idx) = n
                    case (NB_NORTHWEST)
                        northwest(idx) = n
                    case (NB_NORTHEAST)
                        northeast(idx) = n
                    case (NB_SOUTHWEST)
                        southwest(idx) = n
                    case (NB_SOUTHEAST)
                        southeast(idx) = n
                    case default
                        call mpi_exit_on_error("Parcel could not be assigned to a neighbour!")
                end select

                n_sends(nb) = idx + 1
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
                call parcel_deserialize(n_last + n, recvbuf(i:j))
            enddo

            n_last = n_last + recvcount

        end subroutine unpack_parcels


        subroutine allocate_buffers
            if (allocated(north)) then
                return
            endif

            !FIXME

        end subroutine allocate_buffers

end module parcel_mpi
