! =============================================================================
!                       Test parcel laplace correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_laplace_correction
    use unit_test
    use options, only : parcel
    use constants, only : pi, one, zero, two, f14, f32
    use parcel_container, only : n_parcels, parcels, parcel_alloc
    use parcel_bc, only : apply_periodic_bc, apply_reflective_bc
    use parcel_correction, only : apply_laplace             &
                                , lapl_corr_timer           &
                                , init_parcel_correction
    use parcel_interpl, only : vol2grid
    use parcel_ellipse, only : get_ab
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, nz, dx
    use fields, only : volg
    use timer
    implicit none

    double precision :: final_error, init_error, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    double precision :: q, m, delta

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)

    call  parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)

    nx = 32
    nz = 32
    lower  = (/zero, zero/)
    extent = (/one, one/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:nx-1))

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)

    parcel%n_per_cell = 4
    call init_regular_positions

    volg = zero

    parcels%volume = f14 * vcell

    ! b11
    parcels%B(1, :) = get_ab(parcels%volume(1:n_parcels))

    ! b12
    parcels%B(2, :) = zero


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

        call apply_periodic_bc(parcels%position(:, n))
        call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
    enddo

    call vol2grid

    init_error = sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1))

    if (verbose) then
        open(unit=12, file='test_laplace_correction.asc', status='replace')
        write(12,*) '# iteration, average error, max absolute error'
        write(12,*) 0, init_error, maxval(abs(volg(0:nz, 0:nx-1) / vcell - one))
    endif

    if (verbose) then
        open(unit=1235, file='initial_laplace_corr.asc', status='replace')
        write(1235, *) '# x, y volume B11, B22'
        do n = 1, n_parcels
            write(1235, *) parcels%position(:, n), parcels%volume(n), parcels%B(:, n)
        enddo
        close(1235)
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        if (verbose) then
            call vol2grid
            write(12, *) i, sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1)), &
                            maxval(abs(volg(0:nz, 0:nx-1) / vcell - one))
        endif

        if (verbose .and. (i == 2)) then
            open(unit=1236, file='two_iter_laplace_corr.asc', status='replace')
            write(1236, *) '# x, y volume B11, B22'
            do n = 1, n_parcels
                write(1236, *) parcels%position(:, n), parcels%volume(n), parcels%B(:, n)
            enddo
            close(1236)
        endif
    enddo

    if (verbose) then
        close(12)
    endif

    call vol2grid

    final_error = sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1))

    call print_result_dp('Test laplace correction', final_error, init_error)

    deallocate(volg)

end program test_laplace_correction

