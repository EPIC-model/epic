!> This Pencil FFT performs 3D forward and backwards FFTs using pencil decomposition.
!! Borrowed from PMPIC but uses STAFFT (instead of FFTE) for the FFT kernel
!! and this module contains all the data decomposition around this. There is no FFT required
!! in Z, so this performs FFTs in Y and X (in that order forward and reversed backwards).
!! The data decomposition is the complex aspect, there is the concept of forward and backwards
!! transformations. Forward transformations will go from pencil Z to Y to X and the backwards
!! transformations undo these, so go from X to Y to Z.
!! Note that we use quite a lot of buffer space here, this could be cut down if Y=X dimensions
!! so some optimisation on memory could be done there in that case
module fft_pencil
    use mpi_environment
    use mpi_layout
    use constants
    use dimensions
    use stafft
    use sta2dfft
    implicit none

    integer :: Z_INDEX = 1, Y_INDEX = 2, X_INDEX = 3

    !> Describes a specific pencil transposition, from one pencil decomposition to another
    type pencil_layout
        integer :: pencil_size(3)   ! local pencil size (number of grid points)
        integer :: size(3)          ! process decomposition layout
        integer :: coords(3)        ! *this* process location
        integer :: dim              ! pencil dimension
        integer, dimension(:), allocatable :: send_sizes, send_offsets, recv_sizes, recv_offsets
        integer, dimension(:,:), allocatable :: recv_dims, send_dims
    end type pencil_layout

    integer, parameter :: FORWARD=1, BACKWARD=2   !< Transposition directions


    type(communicator) :: fft_y_comm, fft_x_comm !< Communicators for each dimension

    ! Transpositions from one pencil to another
    type(pencil_layout) :: y_from_z_transposition   &
                         , x_from_y_transposition   &
                         , y_from_x_transposition   &
                         , z_from_y_transposition   &
                         , x_from_z_transposition   &
                         , z_from_x_transposition

    ! Temporary buffers used in transposition
    double precision, dimension(:,:,:), contiguous, pointer :: fft_in_y_buffer , fft_in_x_buffer, diff_x_buffer

    logical :: l_initialised = .false.

    integer :: ngrid(3)

    public initialise_pencil_fft, finalise_pencil_fft
contains

    !> Initialises the pencil FFT functionality, this will create the
    !> transposition structures needed
    !! @returns Size of local dimensions in fourier space for this process
    subroutine initialise_pencil_fft(nx, ny, nz)
        integer, intent(in) :: nx, ny, nz
        integer :: x_distinct_sizes(layout%size(I_X)), &
                   y_distinct_sizes(layout%size(I_Y))


        if (l_initialised) then
            return
        endif

        ngrid = (/nz+1, ny, nx/)

        if (layout%l_parallel(I_X) .and. layout%l_parallel(I_Y)) then
            ! Info from https://www.open-mpi.org
            ! Partitions a communicator into subgroups, which form lower-dimensional Cartesian subgrids.
            ! MPI_Cart_sub(comm, remain_dims, newcomm, ierror)
            !   TYPE(MPI_Comm), INTENT(IN) :: comm
            !   LOGICAL, INTENT(IN) :: remain_dims(*)
            !   TYPE(MPI_Comm), INTENT(OUT) :: newcomm
            !   INTEGER, OPTIONAL, INTENT(OUT) :: ierror
            call MPI_Cart_sub(cart%comm, (/.false., .true./), fft_y_comm%comm, cart%err)
            call MPI_Cart_sub(cart%comm, (/.true., .false./), fft_x_comm%comm, cart%err)
            call MPI_Allgather(box%size(I_Y),       &
                               1, MPI_INTEGER,      &
                               y_distinct_sizes,    &
                               1,                   &
                               MPI_INTEGER,         &
                               fft_y_comm%comm,     &
                               fft_y_comm%err)
            call MPI_Allgather(box%size(I_X),       &
                               1,                   &
                               MPI_INTEGER,         &
                               x_distinct_sizes,    &
                               1,                   &
                               MPI_INTEGER,         &
                               fft_x_comm%comm,     &
                               fft_x_comm%err)
        else if (layout%l_parallel(I_Y)) then
            fft_y_comm%comm = cart%comm
            fft_x_comm%comm = MPI_COMM_SELF
            call MPI_Allgather(box%size(I_Y),       &
                               1,                   &
                               MPI_INTEGER,         &
                               y_distinct_sizes,    &
                               1,                   &
                               MPI_INTEGER,         &
                               fft_y_comm%comm,     &
                               fft_y_comm%err)
            x_distinct_sizes = box%size(I_X)
        else if (layout%l_parallel(I_X)) then
            fft_y_comm%comm = MPI_COMM_SELF
            fft_x_comm%comm = cart%comm
            y_distinct_sizes = box%size(I_Y)
            call MPI_Allgather(box%size(I_X),       &
                               1,                   &
                               MPI_INTEGER,         &
                               x_distinct_sizes,    &
                               1,                   &
                               MPI_INTEGER,         &
                               fft_x_comm%comm,     &
                               fft_x_comm%err)
        else
            fft_y_comm%comm = MPI_COMM_SELF
            fft_x_comm%comm = MPI_COMM_SELF
            y_distinct_sizes = box%size(I_Y)
            x_distinct_sizes = box%size(I_X)
        endif

        ! Get communicator sizes:
        call MPI_Comm_size(fft_y_comm%comm, fft_y_comm%size, fft_y_comm%err)
        call MPI_Comm_size(fft_x_comm%comm, fft_x_comm%size, fft_x_comm%err)

        call initialise_transpositions(y_distinct_sizes, x_distinct_sizes)

        call initialise_buffers

        l_initialised = .true.

    end subroutine initialise_pencil_fft

    !> Cleans up allocated buffer memory
    subroutine finalise_pencil_fft
        if ((fft_y_comm%comm .ne. MPI_COMM_SELF) .and. (fft_y_comm%comm .ne. cart%comm)) then
            call MPI_Comm_free(fft_y_comm%comm, cart%err)
        endif

        if ((fft_x_comm%comm .ne. MPI_COMM_SELF) .and. (fft_x_comm%comm .ne. cart%comm)) then
            call MPI_Comm_free(fft_x_comm%comm, cart%err)
        endif

        deallocate(fft_in_y_buffer , fft_in_x_buffer, diff_x_buffer)
    end subroutine finalise_pencil_fft

    !> Initialises memory for the buffers used in the FFT
    subroutine initialise_buffers
        allocate(fft_in_y_buffer(y_from_z_transposition%pencil_size(Y_INDEX),    &
                                 y_from_z_transposition%pencil_size(X_INDEX),    &
                                 y_from_z_transposition%pencil_size(Z_INDEX)))

        allocate(fft_in_x_buffer(x_from_y_transposition%pencil_size(X_INDEX),    &
                                 x_from_y_transposition%pencil_size(Z_INDEX),    &
                                 x_from_y_transposition%pencil_size(Y_INDEX)))

        allocate(diff_x_buffer(x_from_z_transposition%pencil_size(X_INDEX),    &
                               x_from_z_transposition%pencil_size(Z_INDEX),    &
                               x_from_z_transposition%pencil_size(Y_INDEX)))

    end subroutine initialise_buffers

    !> Initialises the pencil transpositions, from a pencil in one dimension to that in another
    !! @param y_distinct_sizes Y sizes per process
    !! @param x_distinct_sizes X sizes per process
    subroutine initialise_transpositions(y_distinct_sizes, x_distinct_sizes)
        integer, dimension(:), intent(in) :: y_distinct_sizes, x_distinct_sizes
        type(pencil_layout)        :: z_pencil

        z_pencil = create_initial_transposition_description(box%size(I_Z))


        ! Transpositions
        y_from_z_transposition = create_transposition(z_pencil, Y_INDEX, y_distinct_sizes, FORWARD)

        x_from_y_transposition = create_transposition(y_from_z_transposition, X_INDEX, &
                                                      x_distinct_sizes, FORWARD)

        y_from_x_transposition=create_transposition(x_from_y_transposition, Y_INDEX, &
                                                    x_distinct_sizes, BACKWARD)

        z_from_y_transposition=create_transposition(y_from_x_transposition, Z_INDEX, &
                                                    y_distinct_sizes, BACKWARD)

    end subroutine initialise_transpositions

    !> Creates a specific pencil transposition description. It is maybe more a decomposition description,
    !! but the main complexity comes from the transposition from existing decomposition to new decomposition so
    !! therefore it is called transposition. The new pencil decomposition depends not only on the dimension to
    !! split on, but also the existing pencil decomposition. The new decomposed dimension (i.e. the existing
    !! pencil dimension) * other local dimensions is used as the sending size, receiving though requires
    !! knowledge about the data size on the source process so others will send
    !! their pencil dimension size to this process.
    !! @param new_pencil_dim The dimension to use as the new pencil decomposition
    !! @param existing_pencil_dim The dimension used in the current decomposition
    !! @param existing_pencil_process_layout Number of processes per dimension for the current decomposition
    !! @param existing_my_location The current processes block location per dimension for the current decomposition
    !! @param existing_pencil_size Pencil size per dimension for the current decomposition
    !! @param process_dim_sizes Sizes of the pencil dimension from other processes that is used
    !!        to calculate receive count
    !! @param direction Whether we are transposing forwards or backwards, backwards is just an inverse
    type(pencil_layout) function create_transposition(existing_transposition,   &
                                                      new_pencil_dim,           &
                                                      process_dim_sizes,        &
                                                      direction)
        type(pencil_layout),   intent(in) :: existing_transposition
        integer, dimension(:), intent(in) :: process_dim_sizes
        integer,               intent(in) :: new_pencil_dim, direction

        create_transposition%size = determine_pencil_process_dimensions(new_pencil_dim,                 &
                                                                        existing_transposition%dim,     &
                                                                        existing_transposition%size)


        create_transposition%coords = determine_my_pencil_location(new_pencil_dim,                  &
                                                                   existing_transposition%dim,      &
                                                                   existing_transposition%coords)

        create_transposition%pencil_size = determine_pencil_size(new_pencil_dim,                &
                                                                 create_transposition%size,     &
                                                                 create_transposition%coords,   &
                                                                 existing_transposition)

        allocate(create_transposition%send_dims(3, create_transposition%size(existing_transposition%dim)), &
                 create_transposition%recv_dims(3, create_transposition%size(existing_transposition%dim)))

        if (direction == FORWARD) then
            call determine_my_process_sizes_per_dim(existing_transposition%dim,         &
                                                    existing_transposition%pencil_size, &
                                                    create_transposition%size,          &
                                                    create_transposition%send_dims)

            call determine_matching_process_dimensions(new_pencil_dim,                      &
                                                       existing_transposition%dim,          &
                                                       process_dim_sizes,                   &
                                                       create_transposition%pencil_size,    &
                                                       create_transposition%size,           &
                                                       create_transposition%recv_dims)
        else
            call determine_my_process_sizes_per_dim(new_pencil_dim,                     &
                                                    create_transposition%pencil_size,   &
                                                    existing_transposition%size,        &
                                                    create_transposition%recv_dims)

            call determine_matching_process_dimensions(existing_transposition%dim,          &
                                                       new_pencil_dim,                      &
                                                       process_dim_sizes,                   &
                                                       existing_transposition%pencil_size,  &
                                                       existing_transposition%size,         &
                                                       create_transposition%send_dims)
        endif

        allocate(create_transposition%send_sizes(size(create_transposition%send_dims, 2)), &
                 create_transposition%send_offsets(size(create_transposition%send_sizes)), &
                 create_transposition%recv_sizes(size(create_transposition%recv_dims, 2)), &
                 create_transposition%recv_offsets(size(create_transposition%recv_sizes)))

        call concatenate_dimension_sizes(create_transposition%send_dims, create_transposition%send_sizes)
        call determine_offsets_from_size(create_transposition%send_sizes, create_transposition%send_offsets)

        call concatenate_dimension_sizes(create_transposition%recv_dims, create_transposition%recv_sizes)
        call determine_offsets_from_size(create_transposition%recv_sizes, create_transposition%recv_offsets)
        create_transposition%dim=new_pencil_dim
    end function create_transposition

    !> Transposes globally to a new pencil decomposition.
    !! This goes from the source dimensions a,b,c to b,c,a (forwards) or c,a,b (backwards).
    !! It requires multiple steps, first the local data is transposed to c,b,a regardless of direction.
    !! then it is communicated via alltoall, each process then assembles its own b,c,a or c,a,b data via
    !! contiguising across blocks as the data layout is nonlinear.
    !! @param transposition_description Description of the transposition
    !! @param source_dims Dimensions of the current pencil that we wish to transpose from, will go from abc to bca
    !! @param comm The MPI communicator associated with the group of processes who will swap data
    !! @param direction Whether this is going forwards or backwards, it makes a difference to the data arrangement
    !! @param source_data Source data (abc)
    !! @param target_data Target data (bca)
    subroutine transpose_to_pencil(transposition_description, source_dims, comm, &
                                   direction, source_data, target_data)
        type(pencil_layout),         intent(in)    :: transposition_description
        integer,                     intent(in)    :: source_dims(3), direction
        type(communicator),          intent(inout) :: comm
        double precision,            intent(in)    :: source_data(:, :, :)
        double precision,            intent(out)   :: target_data(:, :, :)
        double precision, allocatable, save        :: real_temp(:, :, :)
        double precision, allocatable, save        :: real_temp2(:)
        integer                                    :: buf_size

        buf_size = product(transposition_description%pencil_size)

        !$OMP SINGLE
        allocate(real_temp(size(source_data, 3), size(source_data, 2), size(source_data, 1)))
        allocate(real_temp2(buf_size))
        !$OMP END SINGLE

        ! --> realt_temp is x, y, z (c, b, a)
        call rearrange_data_for_sending(real_source=source_data, real_target=real_temp)

        !$OMP SINGLE
        call MPI_Alltoallv(real_temp(1:size(source_data, 3),        &
                                     1:size(source_data, 2),        &
                                     1:size(source_data, 1)),       &
                           transposition_description%send_sizes,    &
                           transposition_description%send_offsets,  &
                           MPI_DOUBLE_PRECISION,                    &
                           real_temp2(1:buf_size),                  &
                           transposition_description%recv_sizes,    &
                           transposition_description%recv_offsets,  &
                           MPI_DOUBLE_PRECISION,                    &
                           comm%comm,                               &
                           comm%err)
        !$OMP END SINGLE

        call contiguise_data(transposition_description,                             &
                             (/source_dims(3), source_dims(2), source_dims(1)/),    &
                             direction,                                             &
                             source_real_buffer=real_temp2,                         &
                             target_real_buffer=target_data)


        !$OMP SINGLE
        deallocate(real_temp, real_temp2)
        !$OMP END SINGLE

    end subroutine transpose_to_pencil

    !> Contiguises from c,b,a to b,c,a (forwards) or c,a,b (backwards) where these are defined by the
    !! source_dims argument. It is not as simple as just swapping the required dimensions, as this is
    !! after the mpi alltoall and each block lies after the previous block running sequentially in a.
    !! @param transposition_description Transposition descriptor
    !! @param source_dims Representation a,b,c of source data, will contiguise to b,a,c
    !! @param direction Whether we wish to contiguise forwards or backwards
    !! @param source_real_buffer Source real data to transform
    !! @param target_real_buffer Target real data which is the result of the operation
    subroutine contiguise_data(transposition_description, source_dims, direction, &
                               source_real_buffer, target_real_buffer)
        integer,                    intent(in)  :: source_dims(3), direction
        type(pencil_layout), intent(in)  :: transposition_description
        double precision,           intent(in)  :: source_real_buffer(:)
        double precision,           intent(out) :: target_real_buffer(:, :, :)
        integer :: number_blocks, i, j, k, n, index_prefix, index_prefix_dim, block_offset, source_index

        number_blocks = size(transposition_description%recv_sizes)
        index_prefix = 0
        block_offset = 0
        index_prefix_dim = merge(2, 1, direction == FORWARD)


        do i = 1, number_blocks
            if (i .ge. 2) then
                index_prefix = index_prefix &
                             + transposition_description%recv_dims(source_dims(index_prefix_dim), i-1)

                block_offset = block_offset + transposition_description%recv_sizes(i-1)
            endif

            !Transformation is either cba -> bca (forward) or cab (backwards)
            !$OMP DO
            do j = 1, transposition_description%recv_dims(source_dims(3), i) ! a
                do k = 1, transposition_description%recv_dims(source_dims(1), i) ! c
                    do n = 1, transposition_description%recv_dims(source_dims(2), i) ! b
                        source_index = block_offset                                                     &
                                    + (j-1) * transposition_description%recv_dims(source_dims(1), i)    &
                                    * transposition_description%recv_dims(source_dims(2), i)            &
                                    + (n-1) * transposition_description%recv_dims(source_dims(1), i)    &
                                    + k

                        if (direction == FORWARD) then
                            target_real_buffer(index_prefix + n, k, j) = source_real_buffer(source_index) ! bca
                        else
                            target_real_buffer(index_prefix + k, j, n) = source_real_buffer(source_index) ! cab
                        endif
                    enddo
                enddo
            enddo
        !$OMP END DO
        enddo

    end subroutine contiguise_data

    !> Rearranges data for sending, transposing a,b,c into c,b,a . This is done as alltoall splits on dimension c
    !! so to go from one pencil to another we assume here that a is the existing pencil as it is contiguous
    !! @param real_source Source data to transpose from
    !! @param real_target Target data to transpose to
    subroutine rearrange_data_for_sending(real_source, real_target)
        double precision, dimension(:, :, :), intent(in)  :: real_source
        double precision, dimension(:, :, :), intent(out) :: real_target
        integer :: i

        !$OMP DO
        do i = 1, size(real_source, 2)
            real_target(:, i, :) = transpose(real_source(:, i, :))
        end do
        !$OMP END DO

    end subroutine rearrange_data_for_sending

    !> Determines the number of elements to on my process per dimension which either need to be sent
    !! to (forwards transformation) or
    !! received from (backwards) each target process (in the row or column)
    !! This depends on the existing pencil decomposition, as effectively we are breaking that contigulity and
    !! decomposing it into n blocks in that dimension now (provided by new_pencil_procs_per_dim)
    !! @param existing_pencil_dim The pencil dimension that we are transforming from
    !! @param existing_pencil_size Existing pencil decomposition sizes per dimension
    !! @param new_pencil_procs_per_dim For the target decomposition the number of processes per dimension
    !! @param global_grid Description of the global grid which we use for sizing information
    subroutine determine_my_process_sizes_per_dim(existing_pencil_dim, existing_pencil_size, &
                                                  new_pencil_procs_per_dim, &
                                                  specific_sizes_per_dim)
        integer, intent(in) :: existing_pencil_dim, existing_pencil_size(:), new_pencil_procs_per_dim(:)
        integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim
        integer :: i, split_size, split_remainder, j, s

        do i=1,3
            if (i == existing_pencil_dim) then
                s = ngrid(i)
                split_size = s / new_pencil_procs_per_dim(i)
                split_remainder = s - split_size * new_pencil_procs_per_dim(i)
                do j = 1, new_pencil_procs_per_dim(existing_pencil_dim)
                    specific_sizes_per_dim(i,j) = merge(split_size+1, split_size, j .le. split_remainder)
                enddo
            else
                specific_sizes_per_dim(i,:) = existing_pencil_size(i)
            endif
        enddo
    end subroutine determine_my_process_sizes_per_dim

    !> Simple helper function to deduce send or receive offsets from the sizes
    !! @param source_sizes Sizes that we are using to build the offsets
    subroutine determine_offsets_from_size(source_sizes, determined_offsets)
        integer, intent(in) :: source_sizes(:)
        integer, dimension(:), intent(inout) :: determined_offsets

        integer :: i

        determined_offsets(1) = 0
        do i = 2, size(source_sizes)
            determined_offsets(i) = determined_offsets(i-1) + source_sizes(i-1)
        enddo
    end subroutine determine_offsets_from_size

    !> Determines the number of processes in each dimension for the target decomposition. This depends heavily
    !! on the existing decomposition, as we basically contiguise our pencil dimension and decompose the existing
    !! pencil dimension. The third dimension remains unchanged
    !! @param new_pencil_dim New pencil dimension
    !! @param existing_pencil_dim Current decomposition pencil dimension
    !! @param existing_pencil_procs Current decomposition process layout
    function determine_pencil_process_dimensions(new_pencil_dim, existing_pencil_dim, existing_pencil_procs)
        integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_pencil_procs(3)
        integer             :: determine_pencil_process_dimensions(3)
        integer             :: i

        do i = 1, 3
            if (i == new_pencil_dim) then
                determine_pencil_process_dimensions(i) = 1
            else if (i == existing_pencil_dim) then
                determine_pencil_process_dimensions(i) = existing_pencil_procs(new_pencil_dim)
            else
                determine_pencil_process_dimensions(i) = existing_pencil_procs(i)
            endif
        enddo
    end function determine_pencil_process_dimensions

    !> Determines my location for each dimension in the new pencil decomposition. I.e. which block I am
    !! operating on
    !! @param new_pencil_dim New pencil decomposition dimension
    !! @param existing_pencil_dim Current pencil dimension
    !! @param existing_locations Location for the current decomposition
    function determine_my_pencil_location(new_pencil_dim, existing_pencil_dim, existing_locations)
        integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_locations(3)
        integer             :: determine_my_pencil_location(3)
        integer             :: i

        do i = 1, 3
            if (i == new_pencil_dim) then
                determine_my_pencil_location(i) = 1
            else if (i == existing_pencil_dim) then
                determine_my_pencil_location(i) = existing_locations(new_pencil_dim)
            else
                determine_my_pencil_location(i) = existing_locations(i)
            endif
        enddo
    end function determine_my_pencil_location

    !> Concatenates sizes in multiple dimensions for each target process (in a row or column) into a product of
    !! that. This represents all the dimension sizes per process
    !! @param dims The sizes, per dimension and per process that we will fold into target process
    subroutine concatenate_dimension_sizes(dims, concatenated_dim_sizes)
        integer, dimension(:,:), intent(in) :: dims
        integer, dimension(:), intent(inout) :: concatenated_dim_sizes
        integer :: i

        do i = 1,size(dims, 2)
            concatenated_dim_sizes(i) = product(dims(:, i))
        enddo
    end subroutine concatenate_dimension_sizes

    !> Determines the sizes per dimension on the matching process either to receive from
    !! (forward transposition) or send to
    !! (backwards transposition) each source process. Not only does this depend on the
    !! my pencil sizes, but it also depends on the amount of data that the source process has to send over
    !! @param new_pencil_dim The dimension for the new pencil decomposition
    !! @param existing_pencil_dim Dimension for the existing pencil decomposition
    !! @param proc_sizes Size of dimension on the source processes (index in array corresponds to source PID)
    !! @param pencil_size My (new) pencil size per dimension
    !! @param pencil_processes_per_dim The process layout per dimension
    subroutine determine_matching_process_dimensions(new_pencil_dim, existing_pencil_dim, proc_sizes, &
                                            pencil_size, pencil_processes_per_dim, specific_sizes_per_dim)
        integer, intent(in) :: new_pencil_dim, existing_pencil_dim, proc_sizes(:), pencil_size(:)
        integer, intent(in) :: pencil_processes_per_dim(:)
        integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim
        integer :: i, j

        do i = 1, pencil_processes_per_dim(existing_pencil_dim)
            do j = 1, 3
                if (j == new_pencil_dim) then
                    specific_sizes_per_dim(j, i) = proc_sizes(i)
                else
                    specific_sizes_per_dim(j, i) = pencil_size(j)
                endif
            enddo
        enddo
    end subroutine determine_matching_process_dimensions

    !> Creates an initial transposition representation of the Z pencil
    !! that MONC is normally decomposed in. This is then
    !! fed into the create transposition procedure which will generate transpositions to other pencils
    type(pencil_layout) function create_initial_transposition_description(nz)
        integer, intent(in) :: nz

        create_initial_transposition_description%dim = Z_INDEX
        create_initial_transposition_description%size(X_INDEX) = layout%size(I_X)
        create_initial_transposition_description%size(Y_INDEX) = layout%size(I_Y)
        create_initial_transposition_description%size(Z_INDEX) = layout%size(I_Z)
        create_initial_transposition_description%coords(X_INDEX) = layout%coords(I_X)
        create_initial_transposition_description%coords(Y_INDEX) = layout%coords(I_Y)
        create_initial_transposition_description%coords(Z_INDEX) = layout%coords(I_Z)
        create_initial_transposition_description%pencil_size(X_INDEX) = box%size(I_X)
        create_initial_transposition_description%pencil_size(Y_INDEX) = box%size(I_Y)
        create_initial_transposition_description%pencil_size(Z_INDEX) = nz
    end function create_initial_transposition_description

    !> Deduces the size of my (local) pencil based upon the new decomposition. This depends heavily on the current
    !! pencil decomposition, the new pencil dimension is the global size, the existing pencil dimension becomes
    !! decomposed based on the number of processes in that dimension. The third dimension remains unchanged.
    !! @param new_pencil_dim Dimension for the new pencil decomposition
    !! @param pencil_process_layout The processes per dimension layout for the new decomposition
    !! @param my_pencil_location My location in the block layout
    !! @param existing_pencil_dim Current decomposition dimension
    !! @param existing_pencil_size Current decomposition sizes
    !! @param global_grid Description of the global grid which we use for sizing information
    function determine_pencil_size(new_pencil_dim, pencil_process_layout, my_pencil_location, &
                                   existing_transposition)

        type(pencil_layout), intent(in) :: existing_transposition
        integer, intent(in) :: new_pencil_dim, pencil_process_layout(3), my_pencil_location(3)
        integer :: determine_pencil_size(3)
        integer :: i, split_size, split_remainder, s

        do i = 1, 3
            if (i == new_pencil_dim) then
                determine_pencil_size(i) = ngrid(new_pencil_dim)
            else if (i == existing_transposition%dim) then
                s = ngrid(i)

                split_size=s/pencil_process_layout(i)
                split_remainder=s - split_size * pencil_process_layout(i)
                determine_pencil_size(i) = merge(split_size+1, split_size, my_pencil_location(i)+1 .le. &
                                                 split_remainder)
            else
                determine_pencil_size(i)=existing_transposition%pencil_size(i)
            endif
        enddo
    end function determine_pencil_size

end module fft_pencil
