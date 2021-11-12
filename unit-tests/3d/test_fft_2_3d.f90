! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module for combined inner
!       dimensions.
! =============================================================================
program test_ft_2_3d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use inversion_utils_mod, only : init_fft &
                                  , xfactors &
                                  , yfactors &
                                  , xtrig    &
                                  , ytrig
    use sta2dfft
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fpp(:, :, :), &
                                     fs1(:, :, :), &
                                     fs2(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z


    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call update_parameters

    allocate(fpp(-1:nz+1, 0:ny-1, 0:nx-1))
    allocate(fs1(0:nz, 0:nx-1, 0:ny-1))
    allocate(fs2(0:nz, 0:nx-1, 0:ny-1))

    call init_fft

    fpp = zero

    ! setup test field
    do i = 0, nx-1
        x = lower(1) + dble(i) * dx(1)
        do j = 0, ny-1
            y = lower(2) + dble(j) * dx(2)
            do k = 0, nz
                z = lower(3) + dble(k) * dx(3)
                fpp(k, j, i) = dcos(four * x) + dsin(y) + dcos(z)
            enddo
        enddo
    enddo

    fs1 = zero
    fs2 = zero


    ! forward FFT
    call combined_inner_dimensions_fft(fpp(0:nz, :, :), fs1)

    call not_combined_inner_dimensions_fft(fpp(0:nz, :, :), fs2)


    ! final check
    error = maxval(dabs(fs1 - fs2))

    call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-15))

    contains
        subroutine combined_inner_dimensions_fft(fp, fs)
            double precision, intent(in)  :: fp(0:nz,ny,nx)  !Physical
            double precision, intent(out) :: fs(0:nz,nx,ny)  !Spectral
            double precision              :: fp_copy(0:nz, ny, nx)
            integer:: kx, iy

            fp_copy = fp

             !Carry out a full x transform first:
            call forfft((nz+1)*ny,nx,fp_copy,xtrig,xfactors)

            !Transpose array:
            do kx=1,nx
                do iy=1,ny
                    fs(:,kx,iy)=fp_copy(:,iy,kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            call forfft((nz+1)*nx,ny,fs,ytrig,yfactors)
        end subroutine combined_inner_dimensions_fft

        subroutine not_combined_inner_dimensions_fft(fp, fs)
            double precision, intent(in)  :: fp(0:nz,ny,nx)  !Physical
            double precision, intent(out) :: fs(0:nz,nx,ny)  !Spectral
            double precision              :: fp_copy(0:nz, ny, nx)
            integer:: kx, iy, iz

            fp_copy = fp

            !Carry out a full x transform first:
            do iz = 0, nz
                call forfft(ny,nx,fp_copy(iz, :, :),xtrig,xfactors)
            enddo

            !Transpose array:
            do kx=1,nx
                do iy=1,ny
                    fs(:,kx,iy)=fp_copy(:,iy,kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            do iz = 0, nz
                call forfft(nx,ny,fs(iz, :, :),ytrig,yfactors)
            enddo
        end subroutine not_combined_inner_dimensions_fft

end program test_ft_2_3d
