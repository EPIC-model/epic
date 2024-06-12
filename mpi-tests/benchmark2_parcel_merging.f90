program benchmark2_parcel_merging
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_layout
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_utils, only : mpi_stop
    use test_utils, only : epic_timer               &
                         , parcel_io_timer          &
                         , register_timer           &
                         , print_timer              &
                         , start_timer              &
                         , stop_timer               &
                         , setup_parcels            &
                         , init_rng
    use parcel_merging
    use parcel_nearest, only : nearest_allreduce_timer, nearest_barrier_timer
    use test_utils, only : merge_nearest_timer      &
                         , merge_tree_resolve_timer
    use parcel_netcdf
    use netcdf_utils
    use netcdf_reader
    use iomanip, only : zfill
    implicit none

    integer        :: allreduce_timer
    character(64)  :: basename
    character(512) :: fname
    integer        :: ncid, n, m, niter
    integer        :: ncells(3), offset, nfiles

    call mpi_env_initialise

    call register_timer('epic', epic_timer)
    call register_timer('MPI allreduce', allreduce_timer)
    call register_timer('parcel I/O', parcel_io_timer)
    call register_timer('parcel merge', merge_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('nearest MPI allreduce', nearest_allreduce_timer)
    call register_timer('nearest MPI barrier', nearest_barrier_timer)

    parcel%lambda_max = 4.0d0
    parcel%min_vratio = 20.0d0
    parcel%size_factor = 1.0d0
    
    call parse_command_line

    fname = trim(basename) // '_' // zfill(offset) // '_parcels.nc'
    call open_netcdf_file(trim(fname), NF90_NOWRITE, ncid)

    call get_netcdf_box(ncid, lower, extent, ncells)
    call close_netcdf_file(ncid)

    nx = ncells(1)
    ny = ncells(2)
    nz = ncells(3)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call parcel_alloc(max_num_parcels)

    call nearest_win_allocate

    call start_timer(epic_timer)

    do n = 0, niter - 1

        m = mod(n, nfiles) + offset

        fname = trim(basename) // '_' // zfill(m) // '_parcels.nc'
        if (world%rank == world%root) then
           print *, "Read:", trim(fname)
        endif
        call read_netcdf_parcels(fname)

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

    call print_timer

    call mpi_env_finalise


    contains
        subroutine parse_command_line
            integer            :: i
            character(len=512) :: arg

            niter = 10
            i = 0
            offset = 0
            nfiles = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--basename') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    basename = trim(arg)
                else if (arg == '--niter') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') niter
                else if (arg == '--offset') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') offset
                else if (arg == '--nfiles') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') nfiles
                else if (arg == '--size_factor') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') parcel%size_factor
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, "./a.out --basename [basename]"
                    endif
                    call mpi_stop
                endif
                i = i+1
            enddo

        end subroutine parse_command_line

end program benchmark2_parcel_merging
