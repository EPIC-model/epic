! =============================================================================
!                        Test gradient correction
!
!         This unit test checks the gradient correction by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 500 relaxation steps.
! =============================================================================
program test_parcel_correction
   use unit_test
    use constants, only : pi, one, zero, f14, f32
    use parcel_container
    use parcel_correction
    use parcel_interpl, only : vol2grid, vol2grid_timer
    use parcel_ellipse, only : get_ab
    use options, only : parcel
    use parameters, only : lower, extent, update_parameters, vcell, dx, nx, nz
    use fields, only : volg
    use timer

    implicit none

    double precision :: final_error, init_error
    integer :: i, j, k, jj, ii

    call  parse_command_line

    call register_timer('vol2grid', vol2grid_timer)
    call register_timer('gradient correction', grad_corr_timer)

    nx = 32
    nz = 32
    lower  = (/-f32, -f32/)
    extent =  (/0.4d0, 0.4d0/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:nx-1))

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)

    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + f14 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + f14 * dx(2) * jj

                    ! add some deviation
                    parcels%position(k, 1) = (one + 1.0e-7) * parcels%position(k, 1)
                    parcels%position(k, 2) = (one + 1.0e-7) * parcels%position(k, 2)
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    volg = zero

    parcel%is_elliptic = .true.

    parcels%volume = f14 * vcell

    ! b11
    parcels%B(:, 1) = get_ab(parcels%volume(1:n_parcels))

    ! b12
    parcels%B(:, 2) = zero

    call vol2grid

    init_error = sum(abs(volg(0:nz, 0:nx-1) - vcell))

    if (verbose) then
        write(*,*) 'Test gradient correction'
        write(*,*) 'iteration, error'
        write(*,*) 0, init_error
    endif

    call init_parcel_correction

    do i = 1, 500
        call vol2grid
        call apply_gradient(volg,1.80d0,0.5d0)
        if (verbose) then
            write(*,*) i, sum(abs(volg(0:nz, 0:nx-1) - vcell))
        endif
    enddo

    final_error = sum(abs(volg(0:nz, 0:nx-1) - vcell))

    call print_result_dp('Test gradient correction', final_error, init_error)

    deallocate(volg)

end program test_parcel_correction
