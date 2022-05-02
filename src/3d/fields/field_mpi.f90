module field_mpi
    use constants, only : zero
    use mpi_layout
    use mpi_communicator
    implicit none

        double precision, allocatable, dimension(:, :, :) :: west_buf,          &
                                                             south_buf,         &
                                                             southwest_buf,     &
                                                             east_halo_buf,     &
                                                             north_halo_buf,    &
                                                             northeast_halo_buf


        double precision, allocatable, dimension(:, :) :: east_buf,             &
                                                          north_buf,            &
                                                          southeast_buf,        &
                                                          northwest_buf,        &
                                                          west_halo_buf,        &
                                                          south_halo_buf,       &
                                                          northwest_halo_buf,   &
                                                          southeast_halo_buf

        double precision, allocatable, dimension(:) :: northeast_buf,       &
                                                       southwest_halo_buf

        integer :: west_size, east_size, north_size, south_size
        integer :: southwest_size, southeast_size, northwest_size, northeast_size

        integer :: west_halo_size, east_halo_size, north_halo_size, south_halo_size
        integer :: southwest_halo_size, southeast_halo_size, northwest_halo_size, northeast_halo_size

        logical :: l_allocated = .false.

        public :: field_halo_reset,          &
                  field_halo_fill,           &
                  field_interior_accumulate, &
                  field_halo_swap

    contains

        subroutine field_halo_reset(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            data(box%lo(3):box%hi(3), box%hlo(2):box%hhi(2),   box%hlo(1))              = zero    ! west halo
            data(box%lo(3):box%hi(3), box%hlo(2):box%hhi(2),   box%hhi(1)-1:box%hhi(1)) = zero    ! east halo
            data(box%lo(3):box%hi(3), box%hlo(2),              box%lo(1):box%hi(1))     = zero    ! south halo
            data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))     = zero    ! north halo

        end subroutine field_halo_reset


        subroutine field_halo_fill(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            call copy_from_interior_to_buffers(data)

            call interior_to_halo_communication

            call copy_from_buffers_to_halo(data, .false.)

        end subroutine field_halo_fill


        subroutine field_interior_accumulate(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            call copy_from_halo_to_buffers(data)

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call copy_from_buffers_to_interior(data, .true.)

        end subroutine field_interior_accumulate


        subroutine field_halo_swap(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_interior_accumulate(data)

            call field_halo_fill(data)

        end subroutine field_halo_swap


        subroutine halo_to_interior_communication
            type(MPI_Request) :: request

            ! send east halo to west buffer
            call MPI_Isend(east_halo_buf, east_halo_size, MPI_DOUBLE, neighbour%east, &
                          HALO_EAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(west_buf, west_size, MPI_DOUBLE, neighbour%west, &
                          HALO_EAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)



            ! send west halo to east buffer
            call MPI_Isend(west_halo_buf, west_halo_size, MPI_DOUBLE, neighbour%west, &
                           HALO_WEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(east_buf, east_size, MPI_DOUBLE, neighbour%east, &
                          HALO_WEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)



            ! send north halo to south buffer
            call MPI_Isend(north_halo_buf, north_halo_size, MPI_DOUBLE, neighbour%north, &
                           HALO_NORTH_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(south_buf, south_size, MPI_DOUBLE, neighbour%south, &
                          HALO_NORTH_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)


            ! send south halo to north buffer
            ! receive north buffer into south halo
            call MPI_Isend(south_halo_buf, south_halo_size, MPI_DOUBLE, neighbour%south, &
                           HALO_SOUTH_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(north_buf, north_size, MPI_DOUBLE, neighbour%north, &
                          HALO_SOUTH_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)



            ! send southwest halo to northeast buffer
            call MPI_Isend(southwest_halo_buf, southwest_halo_size, MPI_DOUBLE, neighbour%southwest, &
                           HALO_SOUTHWEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(northeast_buf, northeast_size, MPI_DOUBLE, neighbour%northeast, &
                          HALO_SOUTHWEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)


            ! send northwest halo to southeast buffer
            call MPI_Isend(northwest_halo_buf, northwest_halo_size, MPI_DOUBLE, neighbour%northwest, &
                           HALO_NORTHWEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(southeast_buf, southeast_size, MPI_DOUBLE, neighbour%southeast, &
                          HALO_NORTHWEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)

            ! send southwest buffer to northeast halo
            call MPI_Isend(northeast_halo_buf, northeast_halo_size, MPI_DOUBLE, neighbour%northeast, &
                           HALO_NORTHEAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(southwest_buf, southwest_size, MPI_DOUBLE, neighbour%southwest, &
                          HALO_NORTHEAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)

            ! send northwest buffer to southeast halo
            call MPI_Isend(southeast_halo_buf, southeast_halo_size, MPI_DOUBLE, neighbour%southeast, &
                           HALO_SOUTHEAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(northwest_buf, northwest_size, MPI_DOUBLE, neighbour%northwest, &
                          HALO_SOUTHEAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)
        end subroutine halo_to_interior_communication


        subroutine interior_to_halo_communication
            type(MPI_Request) :: request

            ! send west buffer to east halo
            call MPI_Isend(west_buf, west_size, MPI_DOUBLE, neighbour%west, &
                           HALO_EAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive west buffer to east halo (left to right)
            call MPI_Recv(east_halo_buf, east_halo_size, MPI_DOUBLE, neighbour%east, &
                          HALO_EAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send east buffer to west halo
            call MPI_Isend(east_buf, east_size, MPI_DOUBLE, neighbour%east, &
                           HALO_WEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive east buffer into west halo (right to left)
            call MPI_Recv(west_halo_buf, west_halo_size, MPI_DOUBLE, neighbour%west, &
                          HALO_WEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)



            ! send south buffer to north halo
            call MPI_Isend(south_buf, south_size, MPI_DOUBLE, neighbour%south, &
                           HALO_NORTH_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive south buffer into north halo
            call MPI_Recv(north_halo_buf, north_halo_size, MPI_DOUBLE, neighbour%north, &
                          HALO_NORTH_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send north buffer to south halo
            call MPI_Isend(north_buf, north_size, MPI_DOUBLE, neighbour%north, &
                           HALO_SOUTH_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive north buffer into south halo
            call MPI_Recv(south_halo_buf, south_halo_size, MPI_DOUBLE, neighbour%south, &
                          HALO_SOUTH_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send northeast buffer to southwest halo
            call MPI_Isend(northeast_buf, northeast_size, MPI_DOUBLE, neighbour%northeast, &
                           HALO_SOUTHWEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive northeast buffer into southwest halo
            call MPI_Recv(southwest_halo_buf, southwest_halo_size, MPI_DOUBLE, neighbour%southwest, &
                          HALO_SOUTHWEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send southeast buffer to northwest halo
            call MPI_Isend(southeast_buf, southeast_size, MPI_DOUBLE, neighbour%southeast, &
                           HALO_NORTHWEST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive southeast buffer into northwest halo
            call MPI_Recv(northwest_halo_buf, northwest_halo_size, MPI_DOUBLE, neighbour%northwest, &
                          HALO_NORTHWEST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send southwest buffer to northeast halo
            call MPI_Isend(southwest_buf, southwest_size, MPI_DOUBLE, neighbour%southwest, &
                           HALO_NORTHEAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive southwest buffer into northeast halo
            call MPI_Recv(northeast_halo_buf, northeast_halo_size, MPI_DOUBLE, neighbour%northeast, &
                           HALO_NORTHEAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)




            ! send northwest buffer to southeast halo
            call MPI_Isend(northwest_buf, northwest_size, MPI_DOUBLE, neighbour%northwest, &
                           HALO_SOUTHEAST_TAG, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            ! receive northwest buffer into southeast halo
            call MPI_Recv(southeast_halo_buf, southeast_halo_size, MPI_DOUBLE, neighbour%southeast, &
                          HALO_SOUTHEAST_TAG, comm_cart, MPI_STATUS_IGNORE, mpi_err)

        end subroutine interior_to_halo_communication


        subroutine allocate_buffers
            integer :: zlen, ylen, xlen

            xlen = box%hi(1)-box%lo(1)+1
            ylen = box%hi(2)-box%lo(2)+1
            zlen = box%hi(3)-box%lo(3)+1

            ! in the context of the sender process
            allocate(west_buf(zlen, ylen, 2))
            allocate(east_buf(zlen, ylen))
            allocate(south_buf(zlen, 2, xlen))
            allocate(north_buf(zlen, xlen))

            allocate(northeast_buf(zlen))
            allocate(southeast_buf(zlen, 2))
            allocate(southwest_buf(zlen, 2, 2))
            allocate(northwest_buf(zlen, 2))

            west_size  = size(west_buf)
            east_size  = size(east_buf)
            south_size = size(south_buf)
            north_size = size(north_buf)

            northeast_size = size(northeast_buf)
            northwest_size = size(northwest_buf)
            southeast_size = size(southeast_buf)
            southwest_size = size(southwest_buf)


            allocate(east_halo_buf(zlen, ylen, 2))
            allocate(west_halo_buf(zlen, ylen))
            allocate(north_halo_buf(zlen, 2, xlen))
            allocate(south_halo_buf(zlen, xlen))

            allocate(southwest_halo_buf(zlen))
            allocate(northwest_halo_buf(zlen, 2))
            allocate(northeast_halo_buf(zlen, 2, 2))
            allocate(southeast_halo_buf(zlen, 2))

            west_halo_size  = size(west_halo_buf)
            east_halo_size  = size(east_halo_buf)
            south_halo_size = size(south_halo_buf)
            north_halo_size = size(north_halo_buf)

            northeast_halo_size = size(northeast_halo_buf)
            northwest_halo_size = size(northwest_halo_buf)
            southeast_halo_size = size(southeast_halo_buf)
            southwest_halo_size = size(southwest_halo_buf)

        end subroutine allocate_buffers

        subroutine copy_from_interior_to_buffers(data)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))

            if (.not. l_allocated) then
                l_allocated = .true.
                call allocate_buffers
            endif

            west_buf  = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf  = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%hi(1))
            south_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%hi(1))
            southeast_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf = data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%lo(1)+1)

        end subroutine copy_from_interior_to_buffers

        subroutine copy_from_halo_to_buffers(data)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))

            if (.not. l_allocated) then
                l_allocated = .true.
                call allocate_buffers
            endif

            east_halo_buf = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1))
            west_halo_buf = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1))
            north_halo_buf = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))
            south_halo_buf = data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1))
            southwest_halo_buf = data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1))
            northwest_halo_buf = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1))
            northeast_halo_buf = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1))
            southeast_halo_buf = data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1))

        end subroutine copy_from_halo_to_buffers

        subroutine copy_from_buffers_to_halo(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add

            if (l_add) then

                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) + east_halo_buf

                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1)) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1)) + west_halo_buf

                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) + north_halo_buf

                data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1)) + south_halo_buf

                data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)) &
                    = data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)) + southwest_halo_buf

                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) &
                    = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) + northwest_halo_buf

                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    + northeast_halo_buf

                data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) + southeast_halo_buf
            else
                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) = east_halo_buf
                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1)) = west_halo_buf
                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) = north_halo_buf
                data(box%lo(3):box%hi(3), box%hlo(2), box%lo(1):box%hi(1)) = south_halo_buf
                data(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)) = southwest_halo_buf
                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) = northwest_halo_buf
                data(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = northeast_halo_buf
                data(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) = southeast_halo_buf
            endif

        end subroutine copy_from_buffers_to_halo

        subroutine copy_from_buffers_to_interior(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add

            if (l_add) then
                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) + west_buf

                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hi(1)) + east_buf

                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) + south_buf

                data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%hi(1)) + north_buf

                data(box%lo(3):box%hi(3), box%hi(2), box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%hi(2), box%hi(1)) + northeast_buf

                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1)) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1)) + southeast_buf

                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) &
                    = data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) + southwest_buf

                data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%lo(1)+1) + northwest_buf

            else
                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1) = west_buf
                data(box%lo(3):box%hi(3), box%lo(2):box%hi(2),   box%hi(1))             = east_buf
                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))   = south_buf
                data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%hi(1))   = north_buf

                data(box%lo(3):box%hi(3), box%hi(2),             box%hi(1))             = northeast_buf
                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1))             = southeast_buf
                data(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = southwest_buf
                data(box%lo(3):box%hi(3), box%hi(2),             box%lo(1):box%lo(1)+1) = northwest_buf
            endif

        end subroutine copy_from_buffers_to_interior

end module field_mpi
