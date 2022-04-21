module mpi_layout
    use mpi_f08
    use mpi_communicator, only : mpi_size, mpi_rank, mpi_err, comm_world, comm_cart
    implicit none

    type box_type
        integer :: lo(3),  hi(3)
        integer :: hlo(3), hhi(3)
    end type box_type

    type neighbour_type
        integer :: west, east
        integer :: south, north
        integer :: southwest, northwest
        integer :: southeast, northeast
    end type neighbour_type

    type(box_type)       :: box
    type(neighbour_type) :: neighbour

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

            if (mpi_size == 1) then
                neighbour%west  = mpi_rank
                neighbour%east = mpi_rank
                neighbour%south = mpi_rank
                neighbour%north = mpi_rank
                box%lo = (/0,  0,  -1  /)
                box%hi = (/nx, ny, nz+1/)
                box%hlo = box%lo
                box%hhi = box%hi
                return
            endif

            ! create slabs, z-direction keeps 1 processor
            dims = (/0, 0/)
            call MPI_Dims_create(mpi_size, 2, dims, mpi_err)

            periods = (/.true., .true./)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
            !   comm_old -- input communicator
            !   ndims    -- number of dimensions of Cartesian grid
            !   dims     -- number of processes in each dimension
            !   periods  -- grid is periodic (true) or not (false) in each dimension
            !   reorder  -- ranking may be reordered (true) or not (false) (logical)
            call MPI_Cart_create(comm_world, 2, dims, periods, .false., comm_cart, mpi_err)

            ! Get MPI rank of corners of local box
            call MPI_Comm_rank(comm_cart, rank, mpi_err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_coords(comm, rank, maxdims, coords, ierror)
            !   comm    -- communicator with Cartesian structure
            !   rank    -- rank of a process within group of comm
            !   maxdims -- length of vector coords in the calling program
            !   coords  -- containing the Cartesian coordinates of the specified process
            call MPI_Cart_coords(comm_cart, rank, 2, coords)

            call set_local_bounds(nx+1, coords(1), dims(1), box%lo(1), box%hi(1))
            call set_local_bounds(ny+1, coords(2), dims(2), box%lo(2), box%hi(2))
            box%lo(3) = 0
            box%hi(3) = nz

            ! box including asymmetric halo
            box%hlo(1:2) = box%lo(1:2) - 1 !- nh - 1
            box%hhi(1:2) = box%hi(1:2) + 2 !nh + 1
            ! we only need 1 halo layer in vertical direction
            box%hlo(3) = -1
            box%hhi(3) = nz + 1

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest)
            !   comm        -- communicator with Cartesian structure
            !   direction   -- coordinate dimension of shift (integer)
            !   disp        -- displacement ( > 0: upward shift, < 0: downward shift) (integer)
            !   rank_source -- rank of source process (integer)
            !   rank_dest   -- rank of destination process (integer)
            call MPI_Cart_shift(comm_cart, 0, 1, neighbour%west,  neighbour%east, mpi_err)
            call MPI_Cart_shift(comm_cart, 1, 1, neighbour%south, neighbour%north, mpi_err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_rank(comm, coords, rank)
            !   comm   -- communicator
            !   coords -- Cartesian coordinates of a process
            !   rank   -- rank of specified process

            ! lower left corner
            call MPI_Cart_rank(comm_cart, (/coords(1)-1, coords(2)-1/), neighbour%southwest, mpi_err)

            ! upper left corner
            call MPI_Cart_rank(comm_cart, (/coords(1)-1, coords(2)+1/), neighbour%northwest, mpi_err)

            ! upper right corner
            call MPI_Cart_rank(comm_cart, (/coords(1)+1, coords(2)+1/), neighbour%northeast, mpi_err)

            ! lower right corner
            call MPI_Cart_rank(comm_cart, (/coords(1)+1, coords(2)-1/), neighbour%southeast, mpi_err)

        end subroutine mpi_layout_init


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
