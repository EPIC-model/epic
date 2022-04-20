module field_mpi
    use mpi_layout
    use mpi_communicator
    implicit none

    integer, parameter :: HALO_WEST_TAG  = 100000,   &
                          HALO_EAST_TAG  = 100001,   &
                          HALO_SOUTH_TAG = 100002,   &
                          HALO_NORTH_TAG = 100003,   &
                          HALO_LOWLE_TAG = 100004,   &
                          HALO_UPPLE_TAG = 100005,   &
                          HALO_UPPRI_TAG = 100006,   &
                          HALO_LOWRI_TAG = 100007

    contains

        subroutine field_halo_fill(data)
            integer :: data(box%hlo(3):box%hhi(3), &
                                     box%hlo(2):box%hhi(2), &
                                     box%hlo(1):box%hhi(1))
            integer :: cnt
            type(MPI_Status) :: status

            ! send to left process -- 2 layers
            cnt = (box%hi(3) - box%lo(3) + 1) * (box%hi(2) - box%lo(2) + 1) * 2

            call MPI_Send(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1), &
                          cnt, MPI_INT, neighbour%west, HALO_WEST_TAG, comm_cart, mpi_err)

            ! receive from right process
            call MPI_Recv(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)), &
                          cnt, MPI_INT, neighbour%east, HALO_WEST_TAG, comm_cart, status, mpi_err)




            ! send to right process -- 1 layer
            cnt = (box%hi(3) - box%lo(3) + 1) * (box%hi(2) - box%lo(2) + 1)

            call MPI_Send(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hi(1)), &
                          cnt, MPI_INT, neighbour%east, HALO_EAST_TAG, comm_cart, mpi_err)

            ! receive from left process
            call MPI_Recv(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1)), &
                          cnt, MPI_INT, neighbour%west, HALO_EAST_TAG, comm_cart, status, mpi_err)



            ! send to lower process -- 2 layers
            cnt = (box%hi(3) - box%lo(3) + 1) * (box%hi(1) - box%lo(1) + 1) * 2

            call MPI_Send(data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)), &
                          cnt, MPI_INT, neighbour%south, HALO_SOUTH_TAG, comm_cart, mpi_err)

            ! receive from upper process
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)), &
                          cnt, MPI_INT, neighbour%north, HALO_SOUTH_TAG, comm_cart, status, mpi_err)





            ! send to upper process -- 1 layer
            cnt = (box%hi(3) - box%lo(3) + 1) * (box%hi(1) - box%lo(1) + 1)

            ! receive from upper process
            call MPI_Send(data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%hi(1)), &
                          cnt, MPI_INT, neighbour%north, HALO_NORTH_TAG, comm_cart, mpi_err)

            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1)), &
                          cnt, MPI_INT, neighbour%south, HALO_NORTH_TAG, comm_cart, status, mpi_err)




            ! send to lower left corner
            cnt = 4

            call MPI_Send(data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1), &
                          cnt, MPI_INT, neighbour%corners(1), HALO_LOWLE_TAG, comm_cart, mpi_err)

            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)), &
                          cnt, MPI_INT, neighbour%corners(3), HALO_LOWLE_TAG, comm_cart, status, mpi_err)




            ! send to upper left corner
            cnt = 2

            call MPI_Send(data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%lo(1)+1), &
                          cnt, MPI_INT, neighbour%corners(2), HALO_UPPLE_TAG, comm_cart, mpi_err)

            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)), &
                          cnt, MPI_INT, neighbour%corners(4), HALO_UPPLE_TAG, comm_cart, status, mpi_err)





            ! send to upper right corner
            cnt = 1

            call MPI_Send(data(box%lo(3):box%hi(3), box%hi(2), box%hi(1)), &
                          cnt, MPI_INT, neighbour%corners(3), HALO_UPPRI_TAG, comm_cart, mpi_err)

            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)), &
                          cnt, MPI_INT, neighbour%corners(1), HALO_UPPRI_TAG, comm_cart, status, mpi_err)



            ! send to lower right corner
            cnt = 2

            call MPI_Send(data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1)), &
                          cnt, MPI_INT, neighbour%corners(4), HALO_LOWRI_TAG, comm_cart, mpi_err)

            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)), &
                          cnt, MPI_INT, neighbour%corners(2), HALO_LOWRI_TAG, comm_cart, status, mpi_err)


        end subroutine field_halo_fill

!         subroutine field_halo_accumulate(data)

!         end subroutine field_halo_accumulate

end module field_mpi
