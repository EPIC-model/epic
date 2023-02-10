module mpi_layout
    use mpi_f08
    use mpi_communicator, only : comm
    implicit none

    type box_type
        integer :: lo(3),  hi(3)
        integer :: hlo(3), hhi(3)
        integer :: size(3)
        integer :: global_ncells(3)
        integer :: ncell
    end type box_type

    type parallel_layout
        logical :: l_parallel(3)
        integer :: size(3)      ! number of processes in each dimension
        integer :: coords(3)    ! Cartesian coordinates of *this* process
    end type parallel_layout

    integer, parameter  :: MPI_NONE      = 0, &
                           MPI_NORTH     = 1, &
                           MPI_SOUTH     = 2, &
                           MPI_WEST      = 3, &
                           MPI_EAST      = 4, &
                           MPI_NORTHWEST = 5, &
                           MPI_NORTHEAST = 6, &
                           MPI_SOUTHWEST = 7, &
                           MPI_SOUTHEAST = 8

    type neighbour_type
        integer          :: rank
        integer          :: lo(2), hi(2)
    end type neighbour_type

    type(box_type)        :: box
    type(neighbour_type)  :: neighbours(8)
    type(parallel_layout) :: layout

    contains

        ! We only distribute x and y.
        ! Each process owns all grid points in z-direction.
        subroutine mpi_layout_init(nx, ny, nz)
            integer, intent(in) :: nx, ny, nz
            integer             :: dims(2)
            integer             :: coords(2)
            integer             :: rank ! we do not reorder the rank numbers, so this is unused!
            logical             :: periods(2)

            ! create slabs, z-direction keeps 1 processor
            dims = (/0, 0/)
            call MPI_Dims_create(comm%size, 2, dims, comm%err)

            layout%size(1) = dims(1)
            layout%size(2) = dims(2)
            layout%size(3) = 1

            periods = (/.true., .true./)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm%cart, ierror)
            !   comm_old -- input communicator
            !   ndims    -- number of dimensions of Cartesian grid
            !   dims     -- number of processes in each dimension
            !   periods  -- grid is periodic (true) or not (false) in each dimension
            !   reorder  -- ranking may be reordered (true) or not (false) (logical)
            call MPI_Cart_create(comm%world, 2, dims, periods, .false., comm%cart, comm%err)

            ! Get MPI rank of corners of local box
            call MPI_Comm_rank(comm%cart, rank, comm%err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_coords(comm, rank, maxdims, coords, ierror)
            !   comm    -- communicator with Cartesian structure
            !   rank    -- rank of a process within group of comm
            !   maxdims -- length of vector coords in the calling program
            !   coords  -- containing the Cartesian coordinates of the specified process
            call MPI_Cart_coords(comm%cart, rank, 2, coords)

            layout%coords(1) = coords(1)
            layout%coords(2) = coords(2)
            layout%coords(3) = 0

            call get_local_bounds(nx, coords(1), dims(1), box%lo(1), box%hi(1))
            call get_local_bounds(ny, coords(2), dims(2), box%lo(2), box%hi(2))
            box%lo(3) = 0
            box%hi(3) = nz

            ! box including asymmetric halo
            box%hlo(1:2) = box%lo(1:2) - 1
            box%hhi(1:2) = box%hi(1:2) + 2
            ! we only need 1 halo layer in vertical direction
            box%hlo(3) = -1
            box%hhi(3) = nz + 1

            layout%l_parallel = (box%hi - box%lo < (/nx-1, ny-1, nz/))
            box%size = box%hi - box%lo + 1
            box%ncell = box%size(1) * box%size(2) * (box%size(3) - 1)

            box%global_ncells = (/nx, ny, nz/)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest)
            !   comm        -- communicator with Cartesian structure
            !   direction   -- coordinate dimension of shift (integer)
            !   disp        -- displacement ( > 0: upward shift, < 0: downward shift) (integer)
            !   rank_source -- rank of source process (integer)
            !   rank_dest   -- rank of destination process (integer)
            call MPI_Cart_shift(comm%cart, 0, 1, neighbours(MPI_WEST)%rank,  neighbours(MPI_EAST)%rank, comm%err)
            call MPI_Cart_shift(comm%cart, 1, 1, neighbours(MPI_SOUTH)%rank, neighbours(MPI_NORTH)%rank, comm%err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_rank(comm, coords, rank)
            !   comm   -- communicator
            !   coords -- Cartesian coordinates of a process
            !   rank   -- rank of specified process

            ! lower left corner
            call MPI_Cart_rank(comm%cart, (/coords(1)-1, coords(2)-1/), neighbours(MPI_SOUTHWEST)%rank, comm%err)

            ! upper left corner
            call MPI_Cart_rank(comm%cart, (/coords(1)-1, coords(2)+1/), neighbours(MPI_NORTHWEST)%rank, comm%err)

            ! upper right corner
            call MPI_Cart_rank(comm%cart, (/coords(1)+1, coords(2)+1/), neighbours(MPI_NORTHEAST)%rank, comm%err)

            ! lower right corner
            call MPI_Cart_rank(comm%cart, (/coords(1)+1, coords(2)-1/), neighbours(MPI_SOUTHEAST)%rank, comm%err)

            ! Obtain limits of neighbours:
            call MPI_Cart_coords(comm%cart, neighbours(MPI_WEST)%rank, 2, coords, comm%err)
            call get_local_bounds(nx, coords(1), dims(1), neighbours(MPI_WEST)%lo(1), neighbours(MPI_WEST)%hi(1))
            neighbours(MPI_WEST)%lo(2) = box%lo(2)
            neighbours(MPI_WEST)%hi(2) = box%hi(2)

            call MPI_Cart_coords(comm%cart, neighbours(MPI_EAST)%rank, 2, coords, comm%err)
            call get_local_bounds(nx, coords(1), dims(1), neighbours(MPI_EAST)%lo(1), neighbours(MPI_EAST)%hi(1))
            neighbours(MPI_EAST)%lo(2) = box%lo(2)
            neighbours(MPI_EAST)%hi(2) = box%hi(2)

            neighbours(MPI_SOUTH)%lo(1) = box%lo(1)
            neighbours(MPI_SOUTH)%hi(1) = box%hi(1)
            call MPI_Cart_coords(comm%cart, neighbours(MPI_SOUTH)%rank, 2, coords, comm%err)
            call get_local_bounds(ny, coords(2), dims(2), neighbours(MPI_SOUTH)%lo(2), neighbours(MPI_SOUTH)%hi(2))

            neighbours(MPI_NORTH)%lo(1) = box%lo(1)
            neighbours(MPI_NORTH)%hi(1) = box%hi(1)
            call MPI_Cart_coords(comm%cart, neighbours(MPI_NORTH)%rank, 2, coords, comm%err)
            call get_local_bounds(ny, coords(2), dims(2), neighbours(MPI_NORTH)%lo(2), neighbours(MPI_NORTH)%hi(2))

            neighbours(MPI_SOUTHWEST)%lo(1) = neighbours(MPI_WEST)%lo(1)
            neighbours(MPI_SOUTHWEST)%hi(1) = neighbours(MPI_WEST)%hi(1)
            neighbours(MPI_SOUTHWEST)%lo(2) = neighbours(MPI_SOUTH)%lo(2)
            neighbours(MPI_SOUTHWEST)%hi(2) = neighbours(MPI_SOUTH)%hi(2)

            neighbours(MPI_NORTHWEST)%lo(1) = neighbours(MPI_WEST)%lo(1)
            neighbours(MPI_NORTHWEST)%hi(1) = neighbours(MPI_WEST)%hi(1)
            neighbours(MPI_NORTHWEST)%lo(2) = neighbours(MPI_NORTH)%lo(2)
            neighbours(MPI_NORTHWEST)%hi(2) = neighbours(MPI_NORTH)%hi(2)

            neighbours(MPI_SOUTHEAST)%lo(1) = neighbours(MPI_EAST)%lo(1)
            neighbours(MPI_SOUTHEAST)%hi(1) = neighbours(MPI_EAST)%hi(1)
            neighbours(MPI_SOUTHEAST)%lo(2) = neighbours(MPI_SOUTH)%lo(2)
            neighbours(MPI_SOUTHEAST)%hi(2) = neighbours(MPI_SOUTH)%hi(2)

            neighbours(MPI_NORTHEAST)%lo(1) = neighbours(MPI_EAST)%lo(1)
            neighbours(MPI_NORTHEAST)%hi(1) = neighbours(MPI_EAST)%hi(1)
            neighbours(MPI_NORTHEAST)%lo(2) = neighbours(MPI_NORTH)%lo(2)
            neighbours(MPI_NORTHEAST)%hi(2) = neighbours(MPI_NORTH)%hi(2)

        end subroutine mpi_layout_init

        pure function is_neighbour(i, j, dir) result(l_inside)
            integer, intent(in) :: i, j, dir
            logical             :: l_inside

            l_inside = ((i >= neighbours(dir)%lo(1)) .and. &
                        (i <= neighbours(dir)%hi(1)) .and. &
                        (j >= neighbours(dir)%lo(2)) .and. &
                        (j <= neighbours(dir)%hi(2)))
        end function is_neighbour

        pure function get_neighbour(i, j) result(nb)
            integer, intent(in) :: i, j
            integer             :: nb, n

            nb = MPI_NONE

            if ((i >= box%lo(1)) .and. &
                (i <= box%hi(1)) .and. &
                (j >= box%lo(2)) .and. &
                (j <= box%hi(2))) then
                return
            endif

            do n = 1, 8
                if (is_neighbour(i, j, n)) then
                    nb = n
                    exit
                endif
            enddo
        end function get_neighbour


        subroutine get_local_bounds(nglobal, coords, dims, first, last)
            integer, intent(in)  :: nglobal, coords, dims
            integer, intent(out) :: first, last
            integer              :: nlocal, remaining

            nlocal = nglobal / dims
            remaining = nglobal - dims * nlocal
            first = nlocal * coords
            if (coords < remaining) then
                nlocal = nlocal + 1
                first = first + coords
                last = last + coords
            else
                first = first + remaining
            endif

            last = first + nlocal - 1

        end subroutine get_local_bounds

end module mpi_layout
