module parcel_mpi
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error
    use mpi_tags
    use fields, only : get_index
    use merge_sort, only : msort
    use parcel_container, only : n_par_attrib       &
                               , parcel_serialize   &
                               , parcel_deserialize &
                               , parcel_delete      &
                               , n_parcels          &
                               , parcels
    implicit none

    private

    ! we have 8 neighbours
    integer :: n_sends(8)
    integer :: n_recvs(8)
    integer :: n_total

    public :: parcel_halo_swap

    contains

        subroutine parcel_halo_swap
            integer :: loc(n_parcels)   !FIXME we need a smart estimate
            integer :: pid(0:n_parcels)   !FIXME we need a smart estimate
            integer :: ii(n_parcels)
            integer :: k

            pid(0) = -1

            ! figure out where parcels should go
            call locate_parcels(loc, pid)

            ! tell your neighbours the number of receiving parcels
            call MPI_Send(n_sends(NB_NORTH), 1, MPI_INT, neighbour%north, NORTH_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_SOUTH), 1, MPI_INT, neighbour%south, NORTH_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_SOUTH), 1, MPI_INT, neighbour%south, SOUTH_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_NORTH), 1, MPI_INT, neighbour%north, SOUTH_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_WEST), 1, MPI_INT, neighbour%west, WEST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_EAST), 1, MPI_INT, neighbour%east, WEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_EAST), 1, MPI_INT, neighbour%east, EAST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_WEST), 1, MPI_INT, neighbour%west, EAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, NORTHWEST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, NORTHWEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_SOUTHEAST), 1, MPI_INT, neighbour%southeast, SOUTHEAST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_NORTHWEST), 1, MPI_INT, neighbour%northwest, SOUTHEAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, NORTHEAST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, NORTHEAST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call MPI_Send(n_sends(NB_SOUTHWEST), 1, MPI_INT, neighbour%southwest, SOUTHWEST_TAG, &
                          comm_cart, mpi_err)

            call MPI_Recv(n_recvs(NB_NORTHEAST), 1, MPI_INT, neighbour%northeast, SOUTHWEST_TAG, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            ! communicate parcels
            call exchange_parcels(pid, n_sends(NB_NORTH), neighbour%north, &
                                       n_recvs(NB_SOUTH), neighbour%south, NORTH_TAG)

            k = n_sends(1)
            call exchange_parcels(pid(k:), n_sends(NB_SOUTH), neighbour%south, &
                                           n_recvs(NB_NORTH), neighbour%north, SOUTH_TAG)

            k = k + n_sends(2)
            call exchange_parcels(pid(k:), n_sends(NB_WEST), neighbour%west, &
                                           n_recvs(NB_EAST), neighbour%east, WEST_TAG)

            k = k + n_sends(3)
            call exchange_parcels(pid(k:), n_sends(NB_EAST), neighbour%east, &
                                           n_recvs(NB_WEST), neighbour%west, EAST_TAG)

            k = k + n_sends(4)
            call exchange_parcels(pid(k:),  &
                                  n_sends(NB_NORTHWEST), neighbour%northwest, &
                                  n_recvs(NB_SOUTHEAST), neighbour%southeast, NORTHWEST_TAG)

            k = k + n_sends(5)
            call exchange_parcels(pid(k:), &
                                  n_sends(NB_NORTHEAST), neighbour%northeast, &
                                  n_recvs(NB_SOUTHWEST), neighbour%southwest, SOUTHEAST_TAG)

            k = k + n_sends(6)
            call exchange_parcels(pid(k:),  &
                                  n_sends(NB_SOUTHWEST), neighbour%southwest, &
                                  n_recvs(NB_NORTHEAST), neighbour%northeast, NORTHEAST_TAG)

            k = k + n_sends(7)
            call exchange_parcels(pid(k:), &
                                  n_sends(NB_SOUTHEAST), neighbour%southeast, &
                                  n_recvs(NB_NORTHWEST), neighbour%northwest, SOUTHWEST_TAG)



            call msort(pid(1:n_total), ii(1:n_total))

            call parcel_delete(pid(0:n_total), n_total)

        end subroutine parcel_halo_swap


        subroutine exchange_parcels(pid, sendcount, dest, recvcount, source, tag)
            integer, intent(in)           :: pid(0:)
            integer, intent(in)           :: sendcount, recvcount
            integer, intent(in)           :: tag
            integer, intent(in)           :: dest, source
            double precision, allocatable :: sendbuf(:), recvbuf(:)
            integer                       :: send_size, recv_size

            send_size = n_par_attrib * sendcount
            recv_size = n_par_attrib * recvcount

            allocate(sendbuf(send_size))
            allocate(recvbuf(recv_size))

            call pack_parcels(pid, sendcount, sendbuf)

            call MPI_Send(sendbuf, send_size, MPI_DOUBLE, dest, tag, comm_cart, mpi_err)

            call MPI_Recv(recvbuf, recv_size, MPI_DOUBLE, source, tag, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)

            call unpack_parcels(recvcount, recvbuf)

            deallocate(sendbuf)
            deallocate(recvbuf)

        end subroutine exchange_parcels

        subroutine locate_parcels(loc, pid)
            integer, intent(out) :: loc(:)
            integer, intent(out) :: pid(0:)
            integer, allocatable :: ii(:), tmp(:)
            integer              :: i, j, k, n, nb, m

            ! reset the number of sends
            n_sends(:) = 0

            m = 1

            do n = 1, n_parcels
                call get_index(parcels%position(:, n), i, j, k)

                nb = get_neighbour(i, j)

                if (nb == NB_NONE) then
                    cycle
                endif

                loc(m) = nb
                pid(m) = n

                m = m + 1

                n_sends(nb) = n_sends(nb) + 1
            enddo

            n_total = sum(n_sends)

            allocate(ii(n_total))
            allocate(tmp(n_total))

            ! sort location in ascending order
            call msort(loc(1:n_total), ii)

            tmp = pid(1:n_total)

            do n = 1, n_total
                pid(n) = tmp(ii(n))
            enddo

            deallocate(ii)
            deallocate(tmp)

        end subroutine locate_parcels


        subroutine pack_parcels(pid, sendcount, sendbuf)
            integer,          intent(in)  :: pid(0:)
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

end module parcel_mpi
