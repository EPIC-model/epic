program benchmark_parcel_merging
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels, vmin
    use parcel_init, only : parcel_default
    use parcel_merging
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_layout
    use mpi_utils, only : mpi_stop
    use test_utils, only : epic_timer               &
                         , merge_timer              &
                         , parcel_io_timer          &
                         , merge_nearest_timer      &
                         , merge_tree_resolve_timer &
                         , init_timer               &
                         , register_timer           &
                         , print_timer              &
                         , start_timer              &
                         , stop_timer
    implicit none

    integer              :: n, k, niter, allreduce_timer, generate_timer, sk
    integer              :: n_orig, n_per_cell
    integer, allocatable :: seed(:)
    double precision     :: rn(12), lam, lam2, abc, a2, b2, c2, theta, phi
    double precision     :: st, ct, sp, cp

    call mpi_env_initialise

    call register_timer('epic', epic_timer)
    call register_timer('parcel merge', merge_timer)
    call register_timer('parcel init', init_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('MPI allreduce', allreduce_timer)
    call register_timer('generate data', generate_timer)

    call parse_command_line

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    parcel%lambda_max = 4.0d0
    parcel%size_factor = 4

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = world%rank
    call random_seed(put=seed)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call start_timer(epic_timer)

    ! make sure each cell has at least some parcels:
    parcel%n_per_cell = 8
    call parcel_default
    n_orig = n_parcels
    parcels%volume(1:n_orig) = 1.2d0 * vmin
    abc = get_abc(parcels%volume(1))
    parcels%B(1, 1:n_orig) = abc ** f23
    parcels%B(4, 1:n_orig) = abc ** f23

    ! account for parcels aready added
    n_per_cell = n_per_cell - 8

    do k = 1, niter

        n_parcels = n_orig + n_per_cell * box%ncell

        ! -------------------------------------------------------------
        ! Set up the parcel configuration:
        call start_timer(generate_timer)

        do n = n_orig+1, n_parcels
            ! rn bewteen 0 and 1
            call random_number(rn)

            ! x = (upper-lower) * r + lower = extent * r + lower
            parcels%position(1, n) = box%extent(1) * rn(1) +  box%lower(1)
            parcels%position(2, n) = box%extent(2) * rn(2) +  box%lower(2)
            parcels%position(3, n) = box%extent(3) * rn(3) +  box%lower(3)

            ! vorticity between -10 and 10: y = 20 * x - 10
            parcels%vorticity(1, n) = 20.0d0 * rn(4) - 10.d0
            parcels%vorticity(2, n) = 20.0d0 * rn(5) - 10.d0
            parcels%vorticity(3, n) = 20.0d0 * rn(6) - 10.d0

            ! buoyancy between -1 and 1: y = 2 * x - 1
            parcels%buoyancy(n) = 2.0d0 * rn(7) - 1.d0

            ! volume between 0.5 * vmin and 1.5 * vmin
            parcels%volume(n) = vmin * rn(8) + f12 * vmin

            ! lam = a / c in [1, 4]
            lam = 3.d0 * rn(9) + 1.0d0

            ! lam2 = a / b
            lam2 = 3.d0 * rn(10) + 1.0d0

            abc = 0.75d0 * parcels%volume(n) / pi

            a2 = (abc * lam * lam2)  ** f23
            b2 = a2 / lam2 ** 2
            c2 = a2 / lam ** 2

            ! theta and phi in [0, 2pi[
            theta = twopi * rn(11)
            phi = twopi * rn(12)

            st = dsin(theta)
            ct = dcos(theta)
            sp = dsin(phi)
            cp = dcos(phi)

            parcels%B(1, n) = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            parcels%B(2, n) = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            parcels%B(3, n) = (a2 - c2) * ct * sp * cp
            parcels%B(4, n) = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            parcels%B(5, n) = (a2 - c2) * st * sp * cp
        enddo

        call parcel_communicate

        call stop_timer(generate_timer)

        n_total_parcels = 0

        call start_timer(allreduce_timer)
        call MPI_Allreduce(n_parcels,       &
                           n_total_parcels, &
                           1,               &
                           MPI_INTEGER,     &
                           MPI_SUM,         &
                           world%comm,      &
                           world%err)

        if (world%rank == world%root) then
            print *, "Number of parcels before merging:", n_total_parcels
        endif

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

        if (world%rank == world%root) then
            print *, "Number of parcels after merging:", n_total_parcels
        endif

        call stop_timer(allreduce_timer)
    enddo

    call stop_timer(epic_timer)

    call nearest_win_deallocate

    call print_timer

    call mpi_env_finalise


    contains
        subroutine parse_command_line
            integer            :: i
            character(len=512) :: arg

            nx = 32
            ny = 32
            nz = 32
            niter = 1
            n_per_cell = 40
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
                else if (arg == '--n_per_cell') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(i6)') n_per_cell
                else if (arg == '--min_vratio') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    read(arg,'(f16.0)') parcel%min_vratio
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, "./a.out --nx [int] --ny [int] --nz [int] ", &
                                 "--niter [int] --n_per_cell [int] --min_vratio [float]"
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
                print *, "n_per_cell", n_per_cell
                print *, "min_vratio", parcel%min_vratio
            endif

        end subroutine parse_command_line

end program benchmark_parcel_merging
