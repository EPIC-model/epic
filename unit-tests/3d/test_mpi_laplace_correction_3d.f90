! =============================================================================
!                       Test 3D parcel laplace correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_mpi_laplace_correction_3d
    use unit_test
    use mpi_environment
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f23, f32, two, four, f12, f18
    use parcel_container, only : n_parcels, parcels, n_total_parcels, parcel_alloc
    use parcel_correction, only : apply_laplace             &
                                , lapl_corr_timer           &
                                , vort_corr_timer           &
                                , init_parcel_correction
    use parcel_interpl, only : vol2grid, halo_swap_timer
    use parcel_ellipsoid, only : get_abc
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, ny, nz, dx
    use fields, only : volg, field_default
    use field_ops
    use parcel_bc
    use parcel_mpi
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_timer
    implicit none

    logical :: passed = .true.
    double precision :: final_error, init_error, max_err, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    integer :: lo(3), hi(3)
    double precision :: q, m, delta

    call mpi_env_initialise

    passed = (world%err == 0)

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = world%rank
    call random_seed(put=seed)

    call parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)
    call register_timer('vorticity correction', vort_corr_timer)
    call register_timer('halo swap', halo_swap_timer)

    nx = 32
    ny = 32
    nz = 32

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    lo = box%lo
    hi = box%hi

    n_parcels = 8*box%ncell
    call parcel_alloc(n_parcels + 1000)

    n_total_parcels = n_parcels
    call mpi_blocking_reduce(n_total_parcels, MPI_SUM, world)

    parcel%n_per_cell = 8
    call init_regular_positions

    ! 0 --> -delta
    ! 1 -->  delta
    ! y = mx + q
    ! -delta = q
    ! delta = m - delta --> m = 2 * delta
    delta = f14 * dx(1)

    q = - delta
    m = two * delta

    ! add some deviation
    do n = 1, n_parcels
        call random_number(val)

        tmp = m * val + q
        parcels%position(1, n) = parcels%position(1, n) + tmp

        call random_number(val)
        tmp = m * val + q
        parcels%position(2, n) = parcels%position(2, n) + tmp

        call random_number(val)
        tmp = m * val + q
        parcels%position(3, n) = parcels%position(3, n) + tmp

        call apply_periodic_bc(parcels%position(:, n))
        call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
    enddo

    call parcel_communicate

    volg = zero

    parcels%volume = f18 * vcell

    parcels%B(:, :) = zero

    ! b11
    parcels%B(1, :) = get_abc(parcels%volume) ** f23

    ! b22
    parcels%B(4, :) = parcels%B(1, :)

    call vol2grid

    volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)

    init_error = get_sum(volg) / (nx * ny * (nz+1))

    max_err = get_abs_max(volg)

    if (l_verbose .and. (world%rank == world%root)) then
        write(*,*) 'test laplace correction'
        write(*,*) 'iteration, average error, max absolute error'
        write(*,*) 0, init_error, max_err
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        if (l_verbose) then
            call vol2grid
            volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)
            final_error = get_sum(volg) / (nx * ny * (nz+1))
            max_err = get_abs_max(volg)
            if (world%rank == world%root) then
                write(*,*) i, final_error, max_err
            endif
        endif
    enddo

    call vol2grid

    volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)

    final_error = get_sum(volg) / (nx * ny * (nz+1))

    passed = (passed .and. (final_error < init_error))

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI laplace correction 3D', passed)
    endif

end program test_mpi_laplace_correction_3d

