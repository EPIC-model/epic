module field_layout
    use mpi_f08
    use mpi_communicator, only : mpi_size, mpi_rank, mpi_err, comm
!     use parameters, only : nx, ny, nz
    implicit none

    type box_type
        integer :: nx, ny, nz
        integer :: lo(3), hi(3)
    end type box_type

    type neighbour_type
        integer :: xlo, xhi
        integer :: ylo, yhi
    end type neighbour_type

    type(box_type)       :: box
    type(neighbour_type) :: neighbour

    private :: set_local_bounds

    contains

        ! We only distribute x and y.
        ! Each process owns all grid points in z-direction.
        subroutine field_layout_init(nx, ny, nz)
            integer, intent(in) :: nx, ny, nz
            integer             :: dims(2)
            type(MPI_Comm)      :: comm_cart
            integer             :: coords(2)
            integer             :: new_rank
            logical             :: periods(2)

            if (mpi_size == 1) then
                neighbour%xlo = mpi_rank
                neighbour%xhi = mpi_rank
                neighbour%ylo = mpi_rank
                neighbour%yhi = mpi_rank
                box%lo = (/0,    0,    0 /)
                box%hi = (/nx-1, ny-1, nz/)
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
            call MPI_Cart_create(comm, 2, dims, periods, .true., comm_cart, mpi_err)

            call MPI_Comm_rank(comm_cart, new_rank, mpi_err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_coords(comm, rank, maxdims, coords, ierror)
            !   comm    -- communicator with Cartesian structure
            !   rank    -- rank of a process within group of comm
            !   maxdims -- length of vector coords in the calling program
            call MPI_Cart_coords(comm_cart, new_rank, 2, coords)

            call set_local_bounds(nx, coords(1), dims(1), box%lo(1), box%hi(1))
            call set_local_bounds(ny, coords(2), dims(2), box%lo(2), box%hi(2))
            box%lo(3) = 0
            box%hi(3) = nz

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest)
            !   comm        -- communicator with Cartesian structure
            !   direction   -- coordinate dimension of shift (integer)
            !   disp        -- displacement ( > 0: upward shift, < 0: downward shift) (integer)
            !   rank_source -- rank of source process (integer)
            !   rank_dest   -- rank of destination process (integer)
            call MPI_Cart_shift(comm_cart, 0, 1, neighbour%xlo, neighbour%xhi, mpi_err)
            call MPI_Cart_shift(comm_cart, 1, 1, neighbour%ylo, neighbour%yhi, mpi_err)

        end subroutine field_layout_init


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

end module field_layout
