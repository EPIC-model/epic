module parcel_nearest_mpi
    use mpi_layout
    use mpi_utils
    use parcel_mpi, only : get_parcel_id_buffer_ptr     &
                         , deallocate_parcel_id_buffers
    use datatypes, only : intlog_pair_t
    use mpi_datatypes, only : MPI_INTEGER_LOGICAL_ARRAY
    implicit none

    private

    logical, allocatable, dimension(:), target :: l_north_buf           &
                                                , l_south_buf           &
                                                , l_west_buf            &
                                                , l_east_buf            &
                                                , l_northwest_buf       &
                                                , l_northeast_buf       &
                                                , l_southwest_buf       &
                                                , l_southeast_buf

    type(intlog_pair_t), allocatable, dimension(:), target :: il_north_buf      &
                                                            , il_south_buf      &
                                                            , il_west_buf       &
                                                            , il_east_buf       &
                                                            , il_northwest_buf  &
                                                            , il_northeast_buf  &
                                                            , il_southwest_buf  &
                                                            , il_southeast_buf

    type :: neighbour_info_t
        integer              :: rank
        integer, allocatable :: iclo(:)       ! store parcel index of remote *ic*
        logical, allocatable :: l_merged(:)
        logical, allocatable :: l_available(:)
        logical, allocatable :: l_leaf(:)

        ! mark all indices that were changed with a *put* operation
        integer, allocatable :: dirty_avail(:)
        integer, allocatable :: dirty_leaf(:)
        integer, allocatable :: dirty_merged(:)

    contains
        procedure :: alloc
        procedure :: dealloc

        procedure, private :: lrealloc
        procedure, private :: irealloc

    end type

    type :: tree_t

        integer :: n_sends(8)
        integer :: n_recvs(8)

        type(neighbour_info_t) :: neighbour_info(8)

        ! Logicals used to determine which mergers are executed
        ! Integers above could be reused for this, but this would
        ! make the algorithm less readable
        logical, allocatable :: l_leaf(:)
        logical, allocatable :: l_available(:)
        logical, allocatable :: l_merged(:)    ! indicates parcels merged in first stage
        integer, allocatable :: index_to_remote(:)  ! entries that must be sent to neighbours

    contains
        procedure, private :: send_from_remote
        procedure, private :: recv_from_remote
        procedure, private :: send_all
        procedure, private :: recv_all

        procedure :: setup

        procedure :: put_avail
        procedure :: put_leaf
        procedure :: put_merged

        procedure :: get_avail
        procedure :: get_leaf
        procedure :: get_merged

        procedure :: sync_avail
        procedure :: sync_leaf
        procedure :: sync_merged

    end type

    public :: tree_t

contains

    subroutine alloc(this, n)
        class(neighbour_info_t), intent(inout) :: this
        integer,                 intent(in)    :: n
        integer                                :: m

        if (.not. allocated(this%l_merged)) then
            allocate(this%l_merged(n))
            allocate(this%l_available(n))
            allocate(this%l_leaf(n))
            allocate(this%iclo(n))
            allocate(this%dirty_avail(n))
            allocate(this%dirty_leaf(n))
            allocate(this%dirty_merged(n))
        else
            m = size(this%l_merged)
            if (n > m) then
                call this%lrealloc(n, m, this%l_merged)
                call this%lrealloc(n, m, this%l_available)
                call this%lrealloc(n, m, this%l_leaf)
                call this%irealloc(n, m, this%iclo)
                call this%irealloc(n, m, this%dirty_avail)
                call this%irealloc(n, m, this%dirty_leaf)
                call this%irealloc(n, m, this%dirty_merged)
            endif
        endif

    end subroutine alloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine lrealloc(this, n, m, l_data)
        class(neighbour_info_t), intent(inout) :: this
        integer,                 intent(in)    :: n
        integer,                 intent(in)    :: m
        logical, allocatable,    intent(inout) :: l_data(:)
        logical, allocatable                   :: tmp(:)

        allocate(tmp(n))
        tmp(1:m) = l_data(1:m)
        call move_alloc(tmp, l_data)

    end subroutine lrealloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine irealloc(this, n, m, l_data)
        class(neighbour_info_t), intent(inout) :: this
        integer,                 intent(in)    :: n
        integer,                 intent(in)    :: m
        integer, allocatable,    intent(inout) :: l_data(:)
        integer, allocatable                   :: tmp(:)

        allocate(tmp(n))
        tmp(1:m) = l_data(1:m)
        call move_alloc(tmp, l_data)

    end subroutine irealloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine dealloc(this)
        class(neighbour_info_t), intent(inout) :: this
        if (allocated(this%l_merged)) then
            deallocate(this%l_merged)
            deallocate(this%l_available)
            deallocate(this%l_leaf)
            deallocate(this%iclo)
            deallocate(this%dirty_avail)
            deallocate(this%dirty_leaf)
            deallocate(this%dirty_merged)
        endif
    end subroutine dealloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            this%l_available(ic) = val
        else
            n = get_neighbour_from_rank(rank)

            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            this%neighbour_info(n)%l_available(m) = val
            this%neighbour_info(n)%dirty_avail(m) = ic
        endif

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            this%l_leaf(ic) = val
        else
            n = get_neighbour_from_rank(rank)
            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            this%neighbour_info(n)%l_leaf(m) = val
            this%neighbour_info(n)%dirty_leaf(m) = ic
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            this%l_merged(ic) = val
        else
            n = get_neighbour_from_rank(rank)

            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            this%neighbour_info(n)%l_merged(m) = val
            this%neighbour_info(n)%dirty_merged(m) = ic
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical                      :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            val = this%l_available(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            val = this%neighbour_info(n)%l_available(m)
        endif

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical                      :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            val = this%l_leaf(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            val = this%neighbour_info(n)%l_leaf(m)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical                      :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            val = this%l_merged(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = findloc(array=this%neighbour_info(n)%iclo, value=ic, dim=1)
            val = this%neighbour_info(n)%l_merged(m)
        endif

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine setup(this, iclo, rclo, n_local_small)
        class(tree_t), intent(inout)   :: this
        integer,       intent(in)      :: iclo(:)
        integer,       intent(in)      :: rclo(:)
        integer,       intent(in)      :: n_local_small
        integer, dimension(:), pointer :: send_buf
        integer, allocatable           :: recv_buf(:)
        type(MPI_Request)              :: requests(8)
        type(MPI_Status)               :: recv_status, send_statuses(8)
        integer                        :: n, m, send_size, recv_size, rc, l

!         call start_timer(nearest_exchange_timer)

        this%n_recvs = 0
        this%n_sends = 0

        !--------------------------------------------------
        ! Figure out how many *this* MPI rank sends to
        ! each neighbour:
        do m = 1, n_local_small
            rc = rclo(m)
            do n = 1, 8
                if (rc == neighbours(n)%rank) then
                    this%n_sends(n) = this%n_sends(n) + 1
!                     exit
                endif
            enddo
        enddo

        ! Send information to neighbouring ranks
        do n = 1, 8
            call MPI_Isend(this%n_sends(n),         &
                           1,                       &
                           MPI_INTEGER,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of parcel_nearest::resolve_tree.")
        enddo

        ! Receive information from neighbouring ranks
        do n = 1, 8
            call MPI_Recv(this%n_recvs(n),          &
                          1,                        &
                          MPI_INTEGER,              &
                          neighbours(n)%rank,       &
                          RECV_NEIGHBOUR_TAG(n),    &
                          cart%comm,                &
                          recv_status,              &
                          cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Recv of parcel_nearest::resolve_tree.")
        enddo

        call MPI_Waitall(8,                 &
                         requests,          &
                         send_statuses,     &
                         cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Waitall of parcel_nearest::resolve_tree.")


        !--------------------------------------------------
        ! Send indices of *ic* to owner remote is pointing
        ! to:

        ! Send information to neighbouring ranks
        do n = 1, 8

            send_size = this%n_sends(n)

            call get_parcel_id_buffer_ptr(n, send_buf)

            allocate(send_buf(send_size))

            l = 1
            do m = 1, n_local_small
                rc = rclo(m)
                if (rc == neighbours(n)%rank) then
                    send_buf(l) = iclo(m)
                    l = l + 1
                endif
            enddo

            call MPI_Isend(send_buf(1:send_size),   &
                           send_size,               &
                           MPI_INTEGER,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of parcel_nearest::resolve_tree.")
        enddo

        ! Receive information from neighbouring ranks
        do n = 1, 8

            recv_size = this%n_recvs(n)

            call this%neighbour_info(n)%alloc(recv_size)

            allocate(recv_buf(recv_size))

            call MPI_Recv(recv_buf(1:recv_size),    &
                          recv_size,                &
                          MPI_INTEGER,              &
                          neighbours(n)%rank,       &
                          RECV_NEIGHBOUR_TAG(n),    &
                          cart%comm,                &
                          recv_status,              &
                          cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Recv of parcel_nearest::resolve_tree.")


            if (recv_size > 0) then
                do m = 1, recv_size
                    this%neighbour_info(n)%iclo(m) = recv_buf(m)
                enddo
            endif

            deallocate(recv_buf)
        enddo

        call MPI_Waitall(8,                 &
                         requests,          &
                         send_statuses,     &
                         cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Waitall of parcel_nearest::setup.")

!         ! Free all send buffers
!         call deallocate_parcel_id_buffers

        print *, "setup done"

!         call stop_timer(nearest_exchange_timer)
    end subroutine setup

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_avail(this)
        class(tree_t), intent(inout) :: this
        type(MPI_Request)            :: requests(8)
        type(MPI_Status)             :: statuses(8)
        integer                      :: n

        !----------------------------------------------------------------------
        ! Send from remote to owning rank and sync data at owning rank
        do n = 1, 8
            call this%send_from_remote(n,                                   &
                                       this%neighbour_info(n)%l_available,  &
                                       this%neighbour_info(n)%dirty_avail,  &
                                       requests(n))
        enddo

        do n = 1, 8
            call this%recv_from_remote(n, this%l_available)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_avail.")

        call free_buffers

        !----------------------------------------------------------------------
        ! Send result from owning rank to remote
        do n = 1, 8
            call this%send_all(n, this%l_available, requests(n))
        enddo

        do n = 1, 8
            call this%recv_all(n, this%neighbour_info(n)%l_available)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_avail.")

        call free_buffers

    end subroutine sync_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_leaf(this)
        class(tree_t), intent(inout) :: this
        type(MPI_Request)            :: requests(8)
        type(MPI_Status)             :: statuses(8)
        integer                      :: n

        !----------------------------------------------------------------------
        ! Send from remote to owning rank and sync data at owning rank
        do n = 1, 8
            call this%send_from_remote(n,                                   &
                                       this%neighbour_info(n)%l_leaf,       &
                                       this%neighbour_info(n)%dirty_leaf,   &
                                       requests(n))
        enddo

        do n = 1, 8
            call this%recv_from_remote(n, this%l_leaf)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_leaf.")

        call free_buffers

        !----------------------------------------------------------------------
        ! Send result from owning rank to remote
        do n = 1, 8
            call this%send_all(n, this%l_leaf, requests(n))
        enddo

        do n = 1, 8
            call this%recv_all(n, this%neighbour_info(n)%l_leaf)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_leaf.")

        call free_buffers

    end subroutine sync_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_merged(this)
        class(tree_t), intent(inout) :: this
        type(MPI_Request)            :: requests(8)
        type(MPI_Status)             :: statuses(8)
        integer                      :: n

        !----------------------------------------------------------------------
        ! Send from remote to owning rank and sync data at owning rank
        do n = 1, 8
            call this%send_from_remote(n,                                   &
                                       this%neighbour_info(n)%l_merged,     &
                                       this%neighbour_info(n)%dirty_merged, &
                                       requests(n))
        enddo

        do n = 1, 8
            call this%recv_from_remote(n, this%l_merged)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_merged.")

        call free_buffers

        print *, "next step"

        !----------------------------------------------------------------------
        ! Send result from owning rank to remote
        do n = 1, 8
            call this%send_all(n, this%l_merged, requests(n))
        enddo

        do n = 1, 8
            call this%recv_all(n, this%neighbour_info(n)%l_merged)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of parcel_mpi::sync_merged.")

        call free_buffers

    end subroutine sync_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine send_from_remote(this, n, l_data, dirty, request)
        class(tree_t),               intent(inout) :: this
        integer,                     intent(in)    :: n
        logical,                     intent(in)    :: l_data(:)
        integer,                     intent(in)    :: dirty(:)
        type(MPI_Request),           intent(inout) :: request
        type(intlog_pair_t), dimension(:), pointer :: send_buf
        integer                                    :: send_size, m, i, l


        call get_int_logical_buffer(n, send_buf)

        ! find chagnged (i.e. 'dirty') values
        send_size = count(dirty > 0)

        print *, "send size:", send_size

        allocate(send_buf(send_size))

        if (send_size > 0) then
            ! pack ic index and logical to send buffer
            l = 1
            do m = 1, size(dirty)
                if (dirty(m) > 0) then
                    i = 1 + (l-1)
                    send_buf(i)%ival = dirty(m)
                    send_buf(i)%lval = l_data(m)
                    l = l + 1
                endif
            enddo
        endif

        call MPI_Isend(send_buf(1:send_size),       &
                       send_size,                   &
                       MPI_INTEGER_LOGICAL_ARRAY,   &
                       neighbours(n)%rank,          &
                       SEND_NEIGHBOUR_TAG(n),       &
                       cart%comm,                   &
                       request,                     &
                       cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Isend of parcel_nearest::send_from_remote.")

    end subroutine send_from_remote

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine recv_from_remote(this, n, l_data)
        class(tree_t),     intent(inout) :: this
        integer,           intent(in)    :: n
        logical,           intent(inout) :: l_data(:)
        type(intlog_pair_t), allocatable :: recv_buf(:)
        type(MPI_Status)                 :: recv_status
        integer                          :: recv_size, m, i, ic


        ! check for incoming messages
        call mpi_check_for_message(neighbours(n)%rank,      &
                                   RECV_NEIGHBOUR_TAG(n),   &
                                   recv_size,               &
                                   cart,                    &
                                   MPI_INTEGER_LOGICAL_ARRAY)


        print *, "recv size:", recv_size

        allocate(recv_buf(recv_size))

        call MPI_Recv(recv_buf(1:recv_size),        &
                      recv_size,                    &
                      MPI_INTEGER_LOGICAL_ARRAY,    &
                      neighbours(n)%rank,           &
                      RECV_NEIGHBOUR_TAG(n),        &
                      cart%comm,                    &
                      recv_status,                  &
                      cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Recv of parcel_nearest::recv_from_remote.")

        if (recv_size > 0) then
            ! unpack ic index and logical to recv buffer
            do m = 1, recv_size
                i = 1 + (m-1)
                ic = recv_buf(i)%ival
                l_data(ic) = recv_buf(i)%lval
            enddo
        endif

        deallocate(recv_buf)

    end subroutine recv_from_remote

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine send_all(this, n, l_data, request)
        class(tree_t),     intent(inout) :: this
        integer,           intent(in)    :: n
        logical,           intent(in)    :: l_data(:)
        type(MPI_Request), intent(inout) :: request
        logical, dimension(:), pointer   :: send_buf
        integer                          :: send_size, m, ic


        call get_logical_buffer(n, send_buf)

        send_size = this%n_sends(n)

        print *, "send size all", send_size, n, SEND_NEIGHBOUR_TAG(n)

        allocate(send_buf(send_size))

        if (send_size > 0) then
            ! pack ic index and logical to send buffer
            do m = 1, send_size
                ic = this%neighbour_info(n)%iclo(m)
                send_buf(m) = l_data(ic)
            enddo
        endif

        call MPI_Isend(send_buf(1:send_size),   &
                       send_size,               &
                       MPI_LOGICAL,             &
                       neighbours(n)%rank,      &
                       SEND_NEIGHBOUR_TAG(n),   &
                       cart%comm,               &
                       request,                 &
                       cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Isend of parcel_nearest::send_all.")

    end subroutine send_all

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine recv_all(this, n, l_data)
        class(tree_t),     intent(inout) :: this
        integer,           intent(in)    :: n
        logical,           intent(inout) :: l_data(:)
        logical, allocatable             :: recv_buf(:)
        type(MPI_Status)                 :: recv_status
        integer                          :: recv_size, m


        ! check for incoming messages
        call mpi_check_for_message(neighbours(n)%rank,      &
                                   RECV_NEIGHBOUR_TAG(n),   &
                                   recv_size,               &
                                   cart,                    &
                                   MPI_LOGICAL)

        print *, "recv size all", recv_size, n, RECV_NEIGHBOUR_TAG(n)

        allocate(recv_buf(recv_size))

        call MPI_Recv(recv_buf(1:recv_size),    &
                      recv_size,                &
                      MPI_LOGICAL,              &
                      neighbours(n)%rank,       &
                      RECV_NEIGHBOUR_TAG(n),    &
                      cart%comm,                &
                      recv_status,              &
                      cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Recv of parcel_nearest::recv_all.")

        if (recv_size > 0) then
            ! unpack ic index and logical to recv buffer
            do m = 1, recv_size
                l_data(m) = recv_buf(m)
            enddo
        endif

        deallocate(recv_buf)

    end subroutine recv_all

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine free_buffers

        if (allocated(l_north_buf)) then
            deallocate(l_north_buf)
            deallocate(l_south_buf)
            deallocate(l_west_buf)
            deallocate(l_east_buf)
            deallocate(l_northwest_buf)
            deallocate(l_northeast_buf)
            deallocate(l_southwest_buf)
            deallocate(l_southeast_buf)
        endif

        if (allocated(il_north_buf)) then
            deallocate(il_north_buf)
            deallocate(il_south_buf)
            deallocate(il_west_buf)
            deallocate(il_east_buf)
            deallocate(il_northwest_buf)
            deallocate(il_northeast_buf)
            deallocate(il_southwest_buf)
            deallocate(il_southeast_buf)
        endif

    end subroutine free_buffers

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine get_logical_buffer(dir, buf_ptr)
        integer,                        intent(in)  :: dir
        logical, dimension(:), pointer, intent(out) :: buf_ptr

        select case (dir)
            case (MPI_NORTH)
                buf_ptr => l_north_buf
            case (MPI_SOUTH)
                buf_ptr => l_south_buf
            case (MPI_WEST)
                buf_ptr => l_west_buf
            case (MPI_EAST)
                buf_ptr => l_east_buf
            case (MPI_NORTHWEST)
                buf_ptr => l_northwest_buf
            case (MPI_NORTHEAST)
                buf_ptr => l_northeast_buf
            case (MPI_SOUTHWEST)
                buf_ptr => l_southwest_buf
            case (MPI_SOUTHEAST)
                buf_ptr => l_southeast_buf
            case default
                call mpi_exit_on_error(&
                    "in parcel_nearest_mpi::get_logical_buffer: No valid direction.")
        end select

    end subroutine get_logical_buffer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine get_int_logical_buffer(dir, buf_ptr)
        integer,                                    intent(in)  :: dir
        type(intlog_pair_t), dimension(:), pointer, intent(out) :: buf_ptr

        select case (dir)
            case (MPI_NORTH)
                buf_ptr => il_north_buf
            case (MPI_SOUTH)
                buf_ptr => il_south_buf
            case (MPI_WEST)
                buf_ptr => il_west_buf
            case (MPI_EAST)
                buf_ptr => il_east_buf
            case (MPI_NORTHWEST)
                buf_ptr => il_northwest_buf
            case (MPI_NORTHEAST)
                buf_ptr => il_northeast_buf
            case (MPI_SOUTHWEST)
                buf_ptr => il_southwest_buf
            case (MPI_SOUTHEAST)
                buf_ptr => il_southeast_buf
            case default
                call mpi_exit_on_error(&
                    "in parcel_nearest_mpi::get_int_logical_buffer: No valid direction.")
        end select

    end subroutine get_int_logical_buffer

end module parcel_nearest_mpi
