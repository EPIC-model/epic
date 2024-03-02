program benchmark_parcel_merging
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels, vmin, dx
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
    integer, allocatable :: seed(:)

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

    call parcel_alloc(max_num_parcels)

    call start_timer(epic_timer)

    do k = 1, niter

        ! -------------------------------------------------------------
        ! Set up the parcel configuration:
        call start_timer(generate_timer)

        call setup_parcels

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

    call parcel_dealloc

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
            parcel%n_per_cell = 40
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
                    read(arg,'(i6)') parcel%n_per_cell
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
                print *, "n_per_cell", parcel%n_per_cell
                print *, "min_vratio", parcel%min_vratio
            endif

        end subroutine parse_command_line

        subroutine setup_parcels
            double precision :: rn(3), lam, lam2, abc, a2, b2, c2, theta, phi
            double precision :: st, ct, sp, cp, corner(3)
            integer          :: ix, iy, iz, m, l

            n_parcels = parcel%n_per_cell * box%ncell

            l = 1
            do iz = 0, nz-1
                do iy = box%lo(2), box%hi(2)
                    do ix = box%lo(1), box%hi(1)
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do m = 1, parcel%n_per_cell
                            ! rn between 0 and 1
                            call random_number(rn)
                            parcels%position(1, l) = corner(1) + dx(1) * rn(1)
                            parcels%position(2, l) = corner(2) + dx(2) * rn(2)
                            parcels%position(3, l) = corner(3) + dx(3) * rn(3)
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_parcels == l - 1) then
                call mpi_exit_on_error("Number of parcels disagree!")
            endif

            do n = 1, n_parcels
                ! rn bewteen 0 and 1
                call random_number(rn)

                ! vorticity between -10 and 10: y = 20 * x - 10
                parcels%vorticity(1, n) = 20.0d0 * rn(1) - 10.d0
                parcels%vorticity(2, n) = 20.0d0 * rn(2) - 10.d0
                parcels%vorticity(3, n) = 20.0d0 * rn(3) - 10.d0

                call random_number(rn)

                ! buoyancy between -1 and 1: y = 2 * x - 1
                parcels%buoyancy(n) = 2.0d0 * rn(1) - 1.d0

                ! volume between 0.5 * vmin and 1.5 * vmin
                parcels%volume(n) = vmin * rn(2) + f12 * vmin

                ! lam = a / c in [1, 4]
                lam = 3.d0 * rn(3) + 1.0d0

                call random_number(rn)

                ! lam2 = a / b
                lam2 = 3.d0 * rn(1) + 1.0d0

                abc = 0.75d0 * parcels%volume(n) / pi

                a2 = (abc * lam * lam2)  ** f23
                b2 = a2 / lam2 ** 2
                c2 = a2 / lam ** 2

                ! theta and phi in [0, 2pi[
                theta = twopi * rn(2)
                phi = twopi * rn(3)

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
        end subroutine setup_parcels

end program benchmark_parcel_merging
