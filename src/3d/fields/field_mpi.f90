module field_mpi
    use mpi_layout
    use mpi_communicator
    implicit none

        double precision, allocatable, dimension(:, :) :: west_buf, east_buf, north_buf, south_buf

        integer :: west_size, east_size, north_size, south_size

        logical :: l_allocated = .false.

        public :: field_halo_fill

    contains

        subroutine field_halo_fill(data)
            integer, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                           box%hlo(2):box%hhi(2), &
                                           box%hlo(1):box%hhi(1))


            call halo_communication(data)

        end subroutine field_halo_fill

!         subroutine field_halo_accumulate(data)

!         end subroutine field_halo_accumulate

        subroutine halo_communication(data)
            integer, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                           box%hlo(2):box%hhi(2), &
                                           box%hlo(1):box%hhi(1))

            integer :: cnt
            type(MPI_Status) :: status

            if (.not. l_allocated) then
                l_allocated = .true.
                call allocate_buffers
            endif

            call copy_to_buffers

            ! send west buffer to east halo
            call MPI_Send(west_buf, west_size, MPI_INT, neighbour%west, HALO_EAST_TAG, comm_cart, mpi_err)

            ! receive west buffer to east halo (left to right)
            call MPI_Recv(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)), &
                          west_size, MPI_INT, neighbour%east, HALO_EAST_TAG, comm_cart, status, mpi_err)




            ! send east buffer to west halo
            call MPI_Send(east_buf, east_size, MPI_INT, neighbour%east, HALO_WEST_TAG, comm_cart, mpi_err)

            ! receive east buffer into west halo (right to left)
            call MPI_Recv(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1)), &
                          east_size, MPI_INT, neighbour%west, HALO_WEST_TAG, comm_cart, status, mpi_err)



            ! send south buffer to north halo
            call MPI_Send(south_buf, south_size, MPI_INT, neighbour%south, HALO_NORTH_TAG, comm_cart, mpi_err)

            ! receive south buffer into north halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)), &
                          south_size, MPI_INT, neighbour%north, HALO_NORTH_TAG, comm_cart, status, mpi_err)




            ! send north buffer to south halo
            call MPI_Send(north_buf, north_size, MPI_INT, neighbour%north, HALO_SOUTH_TAG, comm_cart, mpi_err)

            ! receive north buffer into south halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1)), &
                          north_size, MPI_INT, neighbour%south, HALO_SOUTH_TAG, comm_cart, status, mpi_err)




            ! send northeast buffer to southwest halo
            call MPI_Send(northeast_buf, northeast_size, MPI_INT, neighbour%northeast, &
                          HALO_SOUTHWEST_TAG, comm_cart, mpi_err)

            ! receive northeast buffer into southwest halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)), northeast_size, &
                          MPI_INT, neighbour%southwest, HALO_SOUTHWEST_TAG, comm_cart, status, mpi_err)




            ! send southeast buffer to northwest halo
            call MPI_Send(southeast_buf, &
                          southeast_size, MPI_INT, neighbour%southeast, HALO_NORTHWEST_TAG, comm_cart, mpi_err)

            ! receive southeast buffer into northwest halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)), &
                          cnt, MPI_INT, neighbour%northwest, HALO_NORTHWEST_TAG, comm_cart, status, mpi_err)




            ! send southwest buffer to northeast halo
            call MPI_Send(southwest_buf, southwest_size, &
                          MPI_INT, neighbour%southwest, HALO_NORTHEAST_TAG, comm_cart, mpi_err)

            ! receive southwest buffer into northeast halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)), &
                          southwest_size, MPI_INT, neighbour%northeast, HALO_NORTHEAST_TAG, comm_cart, &
                          status, mpi_err)




            ! send northwest buffer to southeast halo
            call MPI_Send(northwest_buf, northwest_size, &
                          MPI_INT, neighbour%northwest, HALO_SOUTHEAST_TAG, comm_cart, mpi_err)

            ! receive northwest buffer into southeast halo
            call MPI_Recv(data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)), northwest_size, &
                          MPI_INT, neighbour%southeast, HALO_SOUTHEAST_TAG, comm_cart, status, mpi_err)

        end subroutine halo_communication


        subroutine allocate_buffers

            ! in the context of the sender process
            allocate(west_buf( box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1))
            allocate(east_buf( box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%hi(1)))
            allocate(south_buf(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)))
            allocate(north_buf(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%hi(1)))

            allocate(northeast_buf(box%lo(3):box%hi(3), box%hi(2),             box%hi(1)))
            allocate(southeast_buf(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1)))
            allocate(southwest_buf(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1))
            allocate(northwest_buf(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%lo(1)+1))

            west_size  = size(west_buf)
            east_size  = size(east_buf)
            south_size = size(south_buf)
            north_size = size(north_buf)

            northeast_size = size(northeast_buf)
            northwest_size = size(northwest_buf)
            southeast_size = size(southeast_buf)
            southwest_size = size(southwest_buf)

        end subroutine allocate_buffers

        subroutine copy_to_buffers(data)
            integer, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                        box%hlo(2):box%hhi(2), &
                                        box%hlo(1):box%hhi(1))

            west_buf  = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf  = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%hi(1))
            south_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%hi(1))
            southeast_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%lo(1)+1)

        end subroutine copy_to_buffers

end module field_mpi
