! =============================================================================
!                       Test 3D parcel laplace correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_mpi_laplace_correction_3d
    use unit_test
    use mpi_communicator
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f23, f32, two, four, f12
    use parcel_container, only : n_parcels, parcels, parcel_alloc
    use parcel_correction, only : apply_laplace             &
                                , lapl_corr_timer           &
                                , vort_corr_timer           &
                                , init_parcel_correction
    use parcel_interpl, only : vol2grid
    use parcel_ellipsoid, only : get_abc
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, ny, nz
    use fields, only : volg, field_default
    use field_ops
    use mpi_timer
    implicit none

    logical :: passed = .true.
    double precision :: final_error, init_error, max_err, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    double precision, parameter :: dev = 0.005d0
    integer :: lo(3), hi(3)

    call mpi_comm_initialise

    passed = (comm%err == 0)

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = comm%rank
    call random_seed(put=seed)

    call parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)
    call register_timer('vorticity correction', vort_corr_timer)


    nx = 16
    ny = 32
    nz = 64

    lower  = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call field_default

    lo = box%lo
    hi = box%hi

    n_parcels = 27*nz*(hi(1) - lo(1) + 1) * (hi(2) - lo(2) + 1)
    call parcel_alloc(n_parcels + 1000)

    parcel%n_per_cell = 27
    call init_regular_positions

    ! add some deviation
    do n = 1, n_parcels
        do i = 1, 3
            tmp = dev
            call random_number(val)
            if (val < 0.5) then
                tmp = -dev
            endif
            parcels%position(i, n) = parcels%position(i, n) + tmp
        enddo
    enddo

    volg = zero

    parcels%volume = vcell / 27.0d0

    parcels%B(:, :) = zero

    ! b11
    parcels%B(1, :) = get_abc(parcels%volume) ** f23

    ! b22
    parcels%B(4, :) = parcels%B(1, :)

    call vol2grid

    volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)

    init_error = get_sum(volg) / (nx * ny * (nz+1))

    max_err = get_abs_max(volg)

    if (verbose .and. (comm%rank == comm%master)) then
        write(*,*) 'test laplace correction'
        write(*,*) 'iteration, average error, max absolute error'
        write(*,*) 0, init_error, max_err
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        if (verbose) then
            call vol2grid
            volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)
            final_error = get_sum(volg) / (nx * ny * (nz+1))
            max_err = get_abs_max(volg)
            if (comm%rank == comm%master) then
                write(*,*) i, final_error, max_err
            endif
        endif
    enddo

    call vol2grid

    volg(0:nz, lo(2):hi(2), lo(1):hi(1)) = abs(volg(0:nz, lo(2):hi(2), lo(1):hi(1)) / vcell - one)

    final_error = get_sum(volg) / (nx * ny * (nz+1))

    passed = (passed .and. (final_error < init_error))

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI laplace correction 3D', passed)
    endif

end program test_mpi_laplace_correction_3d

