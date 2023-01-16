module mpi_layout
    use mpi_f08
    use mpi_communicator, only : comm
    implicit none

    type box_type
        integer :: lo(3),  hi(3)
        integer :: hlo(3), hhi(3)
        integer :: size(3)
    end type box_type

    type parallel_layout
        logical :: l_parallel(3)
        integer :: size(3)      ! number of processes in each dimension
        integer :: coords(3)    ! Cartesian coordinates of *this* process
    end type parallel_layout

    integer, parameter  :: NB_NONE      = 0, &
                           NB_NORTH     = 1, &
                           NB_SOUTH     = 2, &
                           NB_WEST      = 3, &
                           NB_EAST      = 4, &
                           NB_NORTHWEST = 5, &
                           NB_NORTHEAST = 6, &
                           NB_SOUTHWEST = 7, &
                           NB_SOUTHEAST = 8

    type neighbour_type
        integer :: west, east
        integer :: south, north
        integer :: southwest, northwest
        integer :: southeast, northeast
    end type neighbour_type

    type(box_type)        :: box
    type(neighbour_type)  :: neighbour
    type(parallel_layout) :: layout

    private :: set_local_bounds

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

            call set_local_bounds(nx, coords(1), dims(1), box%lo(1), box%hi(1))
            call set_local_bounds(ny, coords(2), dims(2), box%lo(2), box%hi(2))
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

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest)
            !   comm        -- communicator with Cartesian structure
            !   direction   -- coordinate dimension of shift (integer)
            !   disp        -- displacement ( > 0: upward shift, < 0: downward shift) (integer)
            !   rank_source -- rank of source process (integer)
            !   rank_dest   -- rank of destination process (integer)
            call MPI_Cart_shift(comm%cart, 0, 1, neighbour%west,  neighbour%east, comm%err)
            call MPI_Cart_shift(comm%cart, 1, 1, neighbour%south, neighbour%north, comm%err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_rank(comm, coords, rank)
            !   comm   -- communicator
            !   coords -- Cartesian coordinates of a process
            !   rank   -- rank of specified process

            ! lower left corner
            call MPI_Cart_rank(comm%cart, (/coords(1)-1, coords(2)-1/), neighbour%southwest, comm%err)

            ! upper left corner
            call MPI_Cart_rank(comm%cart, (/coords(1)-1, coords(2)+1/), neighbour%northwest, comm%err)

            ! upper right corner
            call MPI_Cart_rank(comm%cart, (/coords(1)+1, coords(2)+1/), neighbour%northeast, comm%err)

            ! lower right corner
            call MPI_Cart_rank(comm%cart, (/coords(1)+1, coords(2)-1/), neighbour%southeast, comm%err)

        end subroutine mpi_layout_init

        pure function is_north(j) result(l_north)
            integer, intent(in) :: j
            logical             :: l_north

            l_north = (j > box%hi(2))
        end function

        pure function is_south(j) result(l_south)
            integer, intent(in) :: j
            logical             :: l_south

            l_south = (j < box%lo(2))
        end function

        pure function is_west(i) result(l_west)
            integer, intent(in) :: i
            logical             :: l_west

            l_west = (i == box%hlo(1))
        end function

        pure function is_east(i) result(l_east)
            integer, intent(in) :: i
            logical             :: l_east

            l_east = (i > box%hi(1))
        end function

        pure function get_neighbour(i, j) result(nb)
            integer, intent(in) :: i, j
            integer             :: nb

            nb = NB_NONE

            if (is_north(j)) then
                ! check if northwest, north or northeast
                if (is_west(i)) then
                    nb = NB_NORTHWEST
                    return
                else if (is_east(i)) then
                    nb = NB_NORTHEAST
                    return
                endif

                nb = NB_NORTH
                return
            endif

            if (is_south(j)) then
                ! check if southwest, south or southeast
                if (is_west(i)) then
                    nb = NB_SOUTHWEST
                    return
                else if (is_east(i)) then
                    nb = NB_SOUTHEAST
                    return
                endif

                nb = NB_SOUTH
                return
            endif

            ! if none of the above is true, the owner can only be
            ! either neighbour west or east or the rank itself
            if (is_west(i)) then
                nb = NB_WEST
                return
            endif

            if (is_east(i)) then
                nb = NB_EAST
                return
            endif

        end function get_neighbour


        subroutine set_local_bounds(nglobal, coords, dims, first, last)
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

        end subroutine set_local_bounds

end module mpi_layout
