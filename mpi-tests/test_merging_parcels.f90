program test_merging_parcels
    use constants, only : pi, zero, one, two, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_merging
    use parcel_netcdf
    use mpi_environment
    use mpi_layout
    use mpi_utils, only : mpi_stop
    use test_utils, only : epic_timer               &
                         , merge_timer              &
                         , parcel_io_timer          &
                         , merge_nearest_timer      &
                         , merge_tree_resolve_timer &
                         , register_timer           &
                         , print_timer              &
                         , start_timer              &
                         , stop_timer
    implicit none

    integer :: n, niter, allreduce_timer
    logical :: l_write

    call mpi_env_initialise

    call register_timer('epic', epic_timer)
    call register_timer('parcel merge', merge_timer)
    call register_timer('parcel I/O', parcel_io_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('MPI allreduce', allreduce_timer)

    call parse_command_line

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    parcel%lambda_max = 4.0d0
    parcel%size_factor = 4

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call parcel_alloc(max_num_parcels)

    call start_timer(epic_timer)

    do n = 1, niter

        call read_netcdf_parcels('initial_parcels.nc')

        n_total_parcels = 0

        call start_timer(allreduce_timer)
        call MPI_Allreduce(n_parcels,       &
                           n_total_parcels, &
                           1,               &
                           MPI_INTEGER,     &
                           MPI_SUM,         &
                           world%comm,      &
                           world%err)
        call stop_timer(allreduce_timer)

        call parcel_merge

        n_total_parcels = 0
        call start_timer(allreduce_timer)
        call MPI_Allreduce(n_parcels,       &
                           n_total_parcels, &
                           1,               &
                           MPI_INTEGER,     &
                           MPI_SUM,         &
                           world%comm,      &
                           world%err)
        call stop_timer(allreduce_timer)
    enddo


    if (l_write) then
        call create_netcdf_parcel_file('parallel_final', .true., .false.)
        call write_netcdf_parcels(t = 0.0d0)
    endif

    call stop_timer(epic_timer)

    call nearest_win_deallocate

    call print_timer

    call mpi_env_finalise


    contains
        subroutine parse_command_line
            integer            :: i
            character(len=512) :: arg

            l_write = .false.
            nx = 32
            ny = 32
            nz = 32
            niter = 1
            parcel%min_vratio = 40.0d0


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
                else if (arg == '--min_vratio') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') parcel%min_vratio
                else if (arg == '--write-final') then
                    l_write = .true.
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, "./a.out --nx [int] --ny [int] --nz [int] ", &
                                 "--niter [int] --min_vratio [float] (--write-final)"
                    endif
                    call mpi_stop
                endif
                i = i+1
            end do

            if (world%rank == world%root) then
                print *, "nx", nx
                print *, "ny", ny
                print *, "nz", nz
                print *, "niter", niter
                print *, "min_vratio", parcel%min_vratio
            endif

        end subroutine parse_command_line

end program test_merging_parcels
