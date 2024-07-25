program benchmark_parcel_merging
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_init, only : parcel_default
    use parcel_merging
    use parcel_mpi, only : parcel_communicate
    use parcel_nearest, only : nearest_allreduce_timer  &
                             , nearest_barrier_timer    &
                             , nearest_rma_timer
    use mpi_environment
    use mpi_layout
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_utils, only : mpi_stop
    use test_utils, only : epic_timer               &
                         , merge_timer              &
                         , merge_nearest_timer      &
                         , merge_tree_resolve_timer &
                         , register_timer           &
                         , print_timer              &
                         , start_timer              &
                         , stop_timer               &
                         , setup_parcels            &
                         , init_rng
    implicit none

    integer              :: k, niter, allreduce_timer, generate_timer
    double precision     :: lx, ly, lz, xlen, ylen, zlen
    logical              :: l_shuffle, l_variable_nppc
    character(len=1)     :: snum
    integer(kind=int64)  :: buf(9) ! size(n_way_parcel_mergers) = 7; +1 (n_parcel_merges); +1 (n_big_close)

    call mpi_env_initialise

    call register_timer('epic', epic_timer)
    call register_timer('parcel merge', merge_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('MPI allreduce', allreduce_timer)
    call register_timer('generate data', generate_timer)
    call register_timer('nearest MPI allreduce', nearest_allreduce_timer)
    call register_timer('nearest MPI barrier', nearest_barrier_timer)
    call register_timer('nearest MPI RMA', nearest_rma_timer)

    call parse_command_line

    lower  = (/zero, zero, zero/)
    extent = (/lx, ly, lz/)

    parcel%lambda_max = 4.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call parcel_alloc(max_num_parcels)

    call init_rng

    call start_timer(epic_timer)

    do k = 1, niter

        ! -------------------------------------------------------------
        ! Set up the parcel configuration:
        call start_timer(generate_timer)

        call setup_parcels(xlen, ylen, zlen, l_shuffle, l_variable_nppc)

        call stop_timer(generate_timer)

        n_total_parcels = 0

        call start_timer(allreduce_timer)
        call MPI_Allreduce(n_parcels,         &
                           n_total_parcels,   &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)

        call parcel_communicate

        call stop_timer(allreduce_timer)

        if (world%rank == world%root) then
            print *, "Number of parcels before merging:", n_total_parcels
        endif

        call stop_timer(allreduce_timer)

        call parcel_merge

        n_total_parcels = 0
        call start_timer(allreduce_timer)
        call MPI_Allreduce(n_parcels,         &
                           n_total_parcels,   &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)

        if (world%rank == world%root) then
            print *, "Number of parcels after merging:", n_total_parcels
        endif

        call stop_timer(allreduce_timer)
    enddo

    call stop_timer(epic_timer)

    call nearest_win_deallocate

    call parcel_dealloc

    buf(1) = n_parcel_merges
    buf(2) = n_big_close
    buf(3:9) = n_way_parcel_mergers

    call mpi_blocking_reduce(buf, MPI_SUM, world)

    n_parcel_merges = buf(1)
!     n_big_close = buf(2)
    n_way_parcel_mergers = buf(3:9)

    if (world%rank == world%root) then
        print *, "Number of MPI ranks:        ", world%size
        print *, "Total number of merges:     ", n_parcel_merges
        print *, "Number of close big parcels:", buf(2) !n_big_close
        do k = 1, 7
            write(snum, fmt='(I1)')  k+1
            print *, "Number of " // snum // "-way mergers:    ", n_way_parcel_mergers(k)
        enddo
    endif

    call print_timer

    call mpi_env_finalise


    contains
        subroutine parse_command_line
            integer            :: i
            character(len=512) :: arg

            nx = 32
            ny = 32
            nz = 32
            lx = 128.0d0
            ly = 128.0d0
            lz = 128.0d0
            xlen = lx
            ylen = ly
            zlen = lz
            niter = 1
            parcel%n_per_cell = 40
            parcel%min_vratio = 40.0d0
            parcel%size_factor = 1.25d0
            l_shuffle = .false.
            l_variable_nppc = .false.


            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--nx') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') nx
                else if (arg == '--ny') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') ny
                else if (arg == '--nz') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') nz
                else if (arg == '--niter') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') niter
                else if (arg == '--n_per_cell') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') parcel%n_per_cell
                else if (arg == '--min_vratio') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') parcel%min_vratio
                else if (arg == '--size_factor') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') parcel%size_factor
                else if (arg == '--lx') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') lx
                else if (arg == '--ly') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') ly
                else if (arg == '--lz') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') lz
                else if (arg == '--xlen') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') xlen
                else if (arg == '--ylen') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') ylen
                else if (arg == '--zlen') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') zlen
                else if (arg == '--shuffle') then
                    l_shuffle = .true.
                else if (arg == '--variable-nppc') then
                    l_variable_nppc = .true.
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, "./benchmark_parcel_merginga.out ",                &
                                 "--nx [int] --ny [int] --nz [int] ",               &
                                 "--lx [float] --ly [float] --lz [float] ",         &
                                 "--xlen [float] --ylen [float] --zlen [float] ",   &
                                 "--niter [int] --n_per_cell [int] ",               &
                                 "--min_vratio [float] --size_factor [float] ",     &
                                 "--shuffle (optional) --variable-nppc (optional)"
                    endif
                    call mpi_stop
                endif
                i = i+1
            end do

            if (world%rank == world%root) then
                print *, "nx", nx
                print *, "ny", ny
                print *, "nz", nz
                print *, "lx", lx
                print *, "ly", ly
                print *, "lz", lz
                print *, "xlen", xlen
                print *, "ylen", ylen
                print *, "zlen", zlen
                print *, "niter", niter
                print *, "n_per_cell", parcel%n_per_cell
                print *, "min_vratio", parcel%min_vratio
                print *, "size_factor", parcel%size_factor
                print *, "shuffle parcels", l_shuffle
                print *, "variable number of parcels/cell:", l_variable_nppc
            endif

        end subroutine parse_command_line

end program benchmark_parcel_merging
