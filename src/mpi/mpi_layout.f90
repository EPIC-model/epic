module mpi_layout
    use datatypes, only : int64
    use mpi_f08
    use mpi_environment, only : communicator, world
    use mpi_tags
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    type box_type
        integer          :: lo(3),  hi(3)       ! lower and upper grid point
        integer          :: hlo(3), hhi(3)      ! lower and upper halo grid point
        integer          :: size(3)             ! number of grid points excluding halo
        integer          :: halo_size(3)        ! number of grid points including halo
        double precision :: lower(3)            ! origin of *this* subdomain excluding halo
        double precision :: halo_lower(3)       ! origin of *this* subdomain including halo
        double precision :: extent(3)           ! local domain extent excluding halo
        integer          :: global_size(3)    ! global number of cell (excluding halo)
        integer          :: ncell               ! number of cells excluding halo
        integer          :: halo_ncell          ! number of cells including halo
    end type box_type

    type(communicator) :: cart

    type parallel_layout
        logical :: l_parallel(3)
        integer :: size(3)      ! number of processes in each dimension
        integer :: coords(3)    ! Cartesian coordinates of *this* process
    end type parallel_layout

    type neighbour_type
        integer          :: rank
        integer          :: lo(2), hi(2)
    end type neighbour_type

    type(box_type)        :: box
    type(neighbour_type)  :: neighbours(8)
    type(parallel_layout) :: layout

    logical, protected :: l_mpi_layout_initialised = .false.

    double precision, allocatable, dimension(:), target :: mpi_north_buf,       &
                                                           mpi_south_buf,       &
                                                           mpi_west_buf,        &
                                                           mpi_east_buf,        &
                                                           mpi_northwest_buf,   &
                                                           mpi_northeast_buf,   &
                                                           mpi_southwest_buf,   &
                                                           mpi_southeast_buf

    contains

        ! We only distribute x and y.
        ! Each process owns all grid points in z-direction.
        subroutine mpi_layout_init(lower, extent, nx, ny, nz)
            double precision, intent(in) :: lower(3), extent(3)
            integer, intent(in)          :: nx, ny, nz
            integer                      :: dims(2)
            integer                      :: coords(2)
            logical                      :: periods(2)
            double precision             :: dx(3)
            integer(kind=int64)          :: max_size

            if (l_mpi_layout_initialised) then
                return
            endif
            l_mpi_layout_initialised = .true.

            ! create slabs, z-direction keeps 1 processor
            dims = (/0, 0/)
            call MPI_Dims_create(world%size, 2, dims, world%err)

            layout%size(1) = dims(1)
            layout%size(2) = dims(2)
            layout%size(3) = 1

            periods = (/.true., .true./)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_create(comm_old, ndims, dims, periods, reorder, cart%comm, ierror)
            !   comm_old -- input communicator
            !   ndims    -- number of dimensions of Cartesian grid
            !   dims     -- number of processes in each dimension
            !   periods  -- grid is periodic (true) or not (false) in each dimension
            !   reorder  -- ranking may be reordered (true) or not (false) (logical)
            ! Note: This reorder flag is not used although specified by the MPI
            !       specifications. Hence, the rank numbers reorder if set to false
            !       or true. We therefore distinguish between world%rank and cart%rank.
            call MPI_Cart_create(world%comm, 2, dims, periods, .false., cart%comm, world%err)

            ! Get new MPI rank number of *this* process
            call MPI_Comm_rank(cart%comm, cart%rank, cart%err)

            ! Get the number of ranks of the communicator
            call MPI_Comm_size(cart%comm, cart%size, cart%err)

            if (cart%size /= world%size) then
                call mpi_exit_on_error(&
                    "mpi_layout::mpi_layout_init: Wrong communicator size.")
            endif

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_coords(comm, rank, maxdims, coords, ierror)
            !   comm    -- communicator with Cartesian structure
            !   rank    -- rank of a process within group of comm
            !   maxdims -- length of vector coords in the calling program
            !   coords  -- containing the Cartesian coordinates of the specified process
            call MPI_Cart_coords(cart%comm, cart%rank, 2, coords)

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
            box%halo_size = box%hhi - box%hlo + 1
            max_size = box%size(1) * box%size(2) * (box%size(3) - 1)
            if (max_size > huge(box%ncell)) then
                call mpi_exit_on_error(&
                    "Number of local cells larger than an integer can represent. Overflow!")
            endif
            box%ncell = max_size
            box%halo_ncell = box%halo_size(1) * box%halo_size(2) * (box%halo_size(3) - 1)

            box%global_size = (/nx, ny, nz/)
            dx = extent / dble(box%global_size)
            box%lower = lower + dx * dble(box%lo)
            box%halo_lower = box%lower - dx
            box%extent(1) = dx(1) * box%size(1)
            box%extent(2) = dx(2) * box%size(2)
            box%extent(3) = dx(3) * (box%size(3) - 1)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest)
            !   comm        -- communicator with Cartesian structure
            !   direction   -- coordinate dimension of shift (integer)
            !   disp        -- displacement ( > 0: upward shift, < 0: downward shift) (integer)
            !   rank_source -- rank of source process (integer)
            !   rank_dest   -- rank of destination process (integer)
            call MPI_Cart_shift(cart%comm, 0, 1,            &
                                neighbours(MPI_WEST)%rank,  &
                                neighbours(MPI_EAST)%rank,  &
                                cart%err)

            call MPI_Cart_shift(cart%comm, 1, 1,            &
                                neighbours(MPI_SOUTH)%rank, &
                                neighbours(MPI_NORTH)%rank, &
                                cart%err)

            ! Info from https://www.open-mpi.org
            ! MPI_Cart_rank(comm, coords, rank)
            !   comm   -- communicator
            !   coords -- Cartesian coordinates of a process
            !   rank   -- rank of specified process

            ! lower left corner
            call MPI_Cart_rank(cart%comm, (/coords(1)-1, coords(2)-1/), &
                               neighbours(MPI_SOUTHWEST)%rank, cart%err)

            ! upper left corner
            call MPI_Cart_rank(cart%comm, (/coords(1)-1, coords(2)+1/), &
                               neighbours(MPI_NORTHWEST)%rank, cart%err)

            ! upper right corner
            call MPI_Cart_rank(cart%comm, (/coords(1)+1, coords(2)+1/), &
                               neighbours(MPI_NORTHEAST)%rank, cart%err)

            ! lower right corner
            call MPI_Cart_rank(cart%comm, (/coords(1)+1, coords(2)-1/), &
                               neighbours(MPI_SOUTHEAST)%rank, cart%err)

            ! Obtain limits of neighbours:
            call MPI_Cart_coords(cart%comm, neighbours(MPI_WEST)%rank, 2, coords, cart%err)
            call get_local_bounds(nx, coords(1), dims(1),       &
                                  neighbours(MPI_WEST)%lo(1),   &
                                  neighbours(MPI_WEST)%hi(1))
            neighbours(MPI_WEST)%lo(2) = box%lo(2)
            neighbours(MPI_WEST)%hi(2) = box%hi(2)

            call MPI_Cart_coords(cart%comm, neighbours(MPI_EAST)%rank, 2, coords, cart%err)
            call get_local_bounds(nx, coords(1), dims(1),       &
                                  neighbours(MPI_EAST)%lo(1),   &
                                  neighbours(MPI_EAST)%hi(1))
            neighbours(MPI_EAST)%lo(2) = box%lo(2)
            neighbours(MPI_EAST)%hi(2) = box%hi(2)

            neighbours(MPI_SOUTH)%lo(1) = box%lo(1)
            neighbours(MPI_SOUTH)%hi(1) = box%hi(1)
            call MPI_Cart_coords(cart%comm, neighbours(MPI_SOUTH)%rank, 2, coords, cart%err)
            call get_local_bounds(ny, coords(2), dims(2),       &
                                  neighbours(MPI_SOUTH)%lo(2),  &
                                  neighbours(MPI_SOUTH)%hi(2))

            neighbours(MPI_NORTH)%lo(1) = box%lo(1)
            neighbours(MPI_NORTH)%hi(1) = box%hi(1)
            call MPI_Cart_coords(cart%comm, neighbours(MPI_NORTH)%rank, 2, coords, cart%err)
            call get_local_bounds(ny, coords(2), dims(2),       &
                                  neighbours(MPI_NORTH)%lo(2),  &
                                  neighbours(MPI_NORTH)%hi(2))

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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        pure function is_neighbour(i, j, dir) result(l_inside)
            integer, intent(in) :: i, j, dir
            logical             :: l_inside

            l_inside = ((i >= neighbours(dir)%lo(1))  .and. &
                        (i <= neighbours(dir)%hi(1))  .and. &
                        (j >= neighbours(dir)%lo(2))  .and. &
                        (j <= neighbours(dir)%hi(2)))
        end function is_neighbour

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_neighbour(i, j) result(nb)
            integer, intent(in) :: i, j
            integer             :: nb, n
#ifndef NDEBUG
            logical             :: l_found = .false.
#endif
            nb = MPI_NONE

            if ((i >= box%lo(1))  .and. &
                (i <= box%hi(1))  .and. &
                (j >= box%lo(2))  .and. &
                (j <= box%hi(2))) then
                return
            endif

            do n = 1, 8
                if (is_neighbour(i, j, n)) then
                    nb = n
#ifndef NDEBUG
                    l_found = .true.
#endif
                    exit
                endif
            enddo

#ifndef NDEBUG
            if (.not. l_found) then
                call mpi_exit_on_error(&
                    "mpi_layout::get_neighbour: No suitable neighbour found.")
            endif
#endif
        end function get_neighbour

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! As multiple neighbours can have the same rank number,
        ! we cannot uniquely assign a neighbour. We just return
        ! the first neighbour tag we find.
        pure function get_neighbour_from_rank(r) result(nb)
            integer, intent(in) :: r
            integer             :: n, nb

            do n = 1, 8
                if (r == neighbours(n)%rank) then
                    nb = n
                    exit
                endif
            enddo
        end function get_neighbour_from_rank

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_mpi_buffer(dir, buf_ptr)
            integer,                                 intent(in)    :: dir
            double precision, dimension(:), pointer, intent(inout) :: buf_ptr

            select case (dir)
                case (MPI_NORTH)
                    buf_ptr => mpi_north_buf
                case (MPI_SOUTH)
                    buf_ptr => mpi_south_buf
                case (MPI_WEST)
                    buf_ptr => mpi_west_buf
                case (MPI_EAST)
                    buf_ptr => mpi_east_buf
                case (MPI_NORTHWEST)
                    buf_ptr => mpi_northwest_buf
                case (MPI_NORTHEAST)
                    buf_ptr => mpi_northeast_buf
                case (MPI_SOUTHWEST)
                    buf_ptr => mpi_southwest_buf
                case (MPI_SOUTHEAST)
                    buf_ptr => mpi_southeast_buf
                case default
                    call mpi_exit_on_error(&
                        "in mpi_layout::get_mpi_buffer: No valid direction.")
            end select
        end subroutine get_mpi_buffer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_mpi_buffers(buf_sizes)
            integer, intent(in) :: buf_sizes(8)

            allocate(mpi_north_buf(buf_sizes(MPI_NORTH)))
            allocate(mpi_south_buf(buf_sizes(MPI_SOUTH)))
            allocate(mpi_west_buf(buf_sizes(MPI_WEST)))
            allocate(mpi_east_buf(buf_sizes(MPI_EAST)))
            allocate(mpi_northwest_buf(buf_sizes(MPI_NORTHWEST)))
            allocate(mpi_northeast_buf(buf_sizes(MPI_NORTHEAST)))
            allocate(mpi_southwest_buf(buf_sizes(MPI_SOUTHWEST)))
            allocate(mpi_southeast_buf(buf_sizes(MPI_SOUTHEAST)))
        end subroutine allocate_mpi_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine deallocate_mpi_buffers
            if (allocated(mpi_north_buf)) then
                deallocate(mpi_north_buf)
            endif

            if (allocated(mpi_south_buf)) then
                deallocate(mpi_south_buf)
            endif

            if (allocated(mpi_west_buf)) then
                deallocate(mpi_west_buf)
            endif

            if (allocated(mpi_east_buf)) then
                deallocate(mpi_east_buf)
            endif

            if (allocated(mpi_northwest_buf)) then
                deallocate(mpi_northwest_buf)
            endif

            if (allocated(mpi_northeast_buf)) then
                deallocate(mpi_northeast_buf)
            endif

            if (allocated(mpi_southwest_buf)) then
                deallocate(mpi_southwest_buf)
            endif

            if (allocated(mpi_southeast_buf)) then
                deallocate(mpi_southeast_buf)
            endif
        end subroutine deallocate_mpi_buffers

end module mpi_layout
