! =============================================================================
!                       Test 3D parcel laplace correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_laplace_correction_3d
    use unit_test
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f23, f32, two, four, f12
    use parcel_container
    use parcel_correction
    use parcel_interpl, only : vol2grid
    use parcel_ellipsoid, only : get_abc
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, ny, nz
    use fields, only : volg
    use timer
    implicit none

    double precision :: final_error, init_error, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    double precision, parameter :: dev = 0.005d0

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
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

    allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

    n_parcels = 27*nx*ny*nz
    call parcel_alloc(n_parcels)

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

    init_error = sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1))

    if (verbose) then
        write(*,*) 'test laplace correction'
        write(*,*) 'iteration, average error, max absolute error'
        write(*,*) 0, init_error, maxval(abs(volg(0:nz, :, :) / vcell - one))
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        if (verbose) then
            call vol2grid
            write(*,*) i, sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1)), &
                          maxval(abs(volg(0:nz, :, :) / vcell - one))
        endif
    enddo

    call vol2grid

    final_error = sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1))

    call print_result_dp('Test laplace correction 3D', final_error, init_error)

    deallocate(volg)

end program test_laplace_correction_3d

