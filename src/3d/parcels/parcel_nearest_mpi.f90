module parcel_nearest_mpi
    use mpi_layout
    implicit none

    private

    type :: neighbour_info_t
        integer              :: rank
        integer, allocatable :: iclo(:)       ! store parcel incdex of remote *ic*
        logical, allocatable :: l_merged(:)
        logical, allocatable :: l_available(:)
        logical, allocatable :: l_leaf(:)

    contains
        procedure :: alloc
        procedure :: free

        procedure, private :: realloc

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

    end type

    type(tree_t) :: tree


contains

    subroutine alloc(this, n)
        class(neighbour_info_t), intent(inout) :: this
        integer,                 intent(in)    :: n
        integer                                :: m

        if (.not. allocated(this%l_merged)) then
            allocate(this%l_merged)
            allocate(this%l_available)
            allocate(this%l_leaf)
        else
            m = size(this%l_merged)
            if (n > m) then
                call this%realloc(n, m, this%l_merged)
                call this%realloc(n, m, this%l_available)
                call this%realloc(n, m, this%l_leaf)
            endif
        endif

    end subroutine alloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine realloc(this, n, m, l_data)
        class(neighbour_info_t), intent(inout) :: this
        integer,                 intent(in)    :: n
        integer,                 intent(in)    :: m
        logical,                 intent(inout) :: l_data(:)
        logical, allocatable                   :: tmp(:)

        allocate(tmp(n))
        tmp(1:m) = l_data(1:m)
        call move_alloc(tmp, l_data)

    end subroutine realloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine free(this)
        class(neighbour_info_t), intent(inout) :: this
        if (allocated(this%l_merged)) then
            deallocate(this%l_merged)
            deallocate(this%l_available)
            deallocate(this%l_leaf)
        endif
    end subroutine free

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_available(this, rank, ic, val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rank == cart%rank) then
            this%available(ic) = val
        else
            n = get_neighbour_from_rank(rank)

            m = this%neighbour_info(n)%local_index(ic)
            this%neighbour_info(n)%l_available(m) = val
        endif

    end subroutine put_available

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
            m = this%neighbour_info(n)%local_index(ic)
            this%neighbour_info(n)%l_leaf(m) = val
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rc == cart%rank) then
            this%l_merged(ic) = vals
        else
            n = get_neighbour_from_rank(rank)

            m = this%neighbour_info(n)%local_index(ic)
            this%neighbour_info(n)%l_merged(m) = val
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_available(this, rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical                      :: val
        integer                      :: n, m

        if (rc == cart%rank) then
            val = this%l_available(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = this%neighbour_info(n)%local_index(ic)
            val = this%neighbour_info(n)%l_available(m)
        endif

    end function get_available

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical,       intent(in)    :: val
        integer                      :: n, m

        if (rc == cart%rank) then
            val = this%l_leaf(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = this%neighbour_info(n)%local_index(ic)
            val = this%neighbour_info(n)%l_leaf(m)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(rank, ic) result(val)
        class(tree_t), intent(inout) :: this
        integer,       intent(in)    :: rank
        integer,       intent(in)    :: ic
        logical                      :: val
        integer                      :: n, m

        if (rc == cart%rank) then
            val = this%l_merged(ic)
        else
            n = get_neighbour_from_rank(rank)

            m = this%neighbour_info(n)%local_index(ic)
            val = this%neighbour_info(n)%l_merged(m)
        endif

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_avail(this)
        class(tree_t), intent(inout) :: this

        call this%sync(this%l_available)

    end subroutine sync_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_leaf(this)
        class(tree_t), intent(inout) :: this

        call this%sync(this%l_leaf)

    end subroutine sync_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync_merged(this)
        class(tree_t), intent(inout) :: this

        call this%sync(this%l_merged)

    end subroutine sync_merged


    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine setup(this, iclo, rclo, n_local_small)
        class(tree_t), intent(inout)   :: this
        integer,       intent(in)      :: iclo(:)
        integer,       intent(in)      :: rclo(:)
        integer,       intent(in)      :: n_local_small
        integer, dimension(:), pointer :: send_buf
        type(MPI_Request)              :: requests(8)
        type(MPI_Status)               :: recv_status, send_statuses(8)
        integer                        :: n, m, send_size, recv_size

        call start_timer(nearest_exchange_timer)

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
                    exit
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

            call get_parcel_id_ptr(n, send_buf)

            allocate(send_buf(send_size))

            do m = 1, send_size
                send_buf(m) = this%neighbour_info(n)%iclo(m)
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

            this%neighbour_info(n)%alloc(recv_size)

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

        enddo

        call MPI_Waitall(8,                 &
                         requests,          &
                         send_statuses,     &
                         cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Waitall of parcel_nearest::setup.")

        ! Free all send buffers
        call deallocate_parcel_id_buffers

        call stop_timer(nearest_exchange_timer)
    end subroutine setup

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine sync(this, l_remote, l_owner)
        class(tree_t), intent(inout) :: this
        logical,       intent(in)    :: l_remote(:)
        logical,       intent(inout) :: l_owner(:)
        integer                      :: n, n_entries

        ! 1. send from remote to owning rank
        ! 2. sync data at owning rank
        ! 3. send result from owning rank to remote

        ! we only send logical + index
        n_entries = 2
        call allocate_parcel_buffers(n_entries)


        call this%send_from_remote_to_owner(l_remote, l_owner)

        call this%send_from_owner_to_remote(l_remote, l_owner)

    end subroutine sync


    subroutine send_from_remote_to_owner(this, l_remote, l_owner)
        class(tree_t), intent(inout)            :: this
        logical,       intent(in)               :: l_remote(:)
        logical,       intent(inout)            :: l_owner(:)
        integer,          dimension(:), pointer :: send_ptr
        double precision, dimension(:), pointer :: send_buf
        double precision, allocatable           :: recv_buf(:)
        type(MPI_Request)                       :: requests(8)
        type(MPI_Status)                        :: recv_status, send_statuses(8)
        integer                                 :: n, send_size, n_entries, m, i


        do n = 1, 8

            call get_parcel_buffer_ptr(n, send_ptr, send_buf)

            send_size = n_sends(n) * n_entries

            if (n_sends(n) > 0) then
                ! pack ic index and logical to send buffer
                do m = 1, n_sends(n)
                    i = 1 + (m-1) * n_entries
                    ic = this%neighbour_info(n)%iclo(m)

                    send_buf(i) = dble(ic)
                    if (l_data(m)) then
                        send_buf(i+1) = 1.0d0
                    else
                        send_buf(i+1) = 0.0d0
                    endif
                enddo
            endif

            call MPI_Isend(send_buf(1:send_size),   &
                           send_size,               &
                           MPI_DOUBLE_PRECISION,    &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of parcel_nearest::exchange_tree_data.")
        enddo


        do n = 1, 8

            ! check for incoming messages
            call mpi_check_for_message(neighbours(n)%rank,      &
                                       RECV_NEIGHBOUR_TAG(n),   &
                                       recv_size,               &
                                       cart)

            allocate(recv_buf(recv_size))

            call MPI_Recv(recv_buf(1:recv_size),    &
                          recv_size,                &
                          MPI_DOUBLE_PRECISION,     &
                          neighbours(n)%rank,       &
                          RECV_NEIGHBOUR_TAG(n),    &
                          cart%comm,                &
                          recv_status,              &
                          cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Recv of parcel_nearest::exchange_tree_data.")

            if (mod(recv_size, n_entries) /= 0) then
                call mpi_exit_on_error(&
                    "parcel_nearest::exchange_tree_data: Receiving wrong count.")
            endif

            recv_count = recv_size / n_entries

            if (recv_count > 0) then
                ! unpack ic index and logical to recv buffer
                do m = 1, recv_count
                    i = 1 + (m-1) * n_entries
                    ic = recv_buf(i)
                    l_owner(ic) = (recv_buf(i+1) > 0.0d0)
                enddo
            endif

            deallocate(recv_buf)
        enddo

    end subroutine send_from_remote_to_owner

end module parcel_nearest_mpi
