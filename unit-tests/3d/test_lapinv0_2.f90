! =============================================================================
!                       Test subroutine lapinv0
!
!  This unit test checks the subroutine lapinv0.
! =============================================================================
program test_lapinv0_2
    use unit_test
    use options, only : verbose
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent, ncell
    use inversion_utils, only : init_fft, fftxyp2s, fftxys2p
    use inversion_mod, only : lapinv0
    implicit none

    !Gridded x component of vorticity (note array order in 3D):
    double precision, allocatable :: vorx(:, :, :)

    !Corresponding potential and exact form:
    double precision, allocatable :: potx(:, :, :), potxe(:, :, :)

    !Spectral work array:
    double precision, allocatable :: ds(:, :, :)

    double precision :: akx, aky, akz  !k, l, m
    double precision :: sikx, siky, sikz
    double precision :: vmult, erms, emax
    integer :: iz, iy, ix

    call parse_command_line

    nx = 32
    ny = 64
    nz = 128
    lower  = (/zero, zero, zero/)
    extent =  (/pi, twopi, two * twopi/)

    call update_parameters

    allocate(vorx(0:nz, ny, nx))
    allocate(potx(0:nz, ny, nx))
    allocate(potxe(0:nz, ny, nx))
    allocate(ds(0:nz, nx, ny))

    !-------------------------------------------------------
    ! Initialise
    call init_fft

    ! Define exact potential and vorticity:
    akx = twopi / extent(1)
    aky = twopi / extent(2)
    akz =    pi / extent(3)
    vmult = -(akx ** 2 + aky ** 2 + akz ** 2)

    do ix = 1, nx
        sikx = dsin(akx * (lower(1) + dx(1) * dble(ix-1)))
        do iy = 1, ny
            siky = dsin(aky * (lower(2) + dx(2) *dble(iy-1)))
            do iz = 0, nz
                sikz = dsin(akz * dx(3) * dble(iz))
                potxe(iz, iy, ix) = sikx * siky * sikz
                vorx(iz, iy, ix) = vmult * potxe(iz, iy, ix)
            enddo
        enddo
    enddo

    !-------------------------------------------------------
    ! Invert vorticity to obtain potential:
    call fftxyp2s(vorx, ds)
    call lapinv0(ds)
    call fftxys2p(ds, potx)

    !-------------------------------------------------------
    ! Check error:
    potx = potx - potxe
    emax = maxval(abs(potx))
    erms = sqrt((f12 * sum(potx(0, :, :) ** 2 + potx(nz,:,:) ** 2) &
                     + sum(potx(1:nz-1, :, :) ** 2) ) / dble(ncell))

    if (verbose) then
        write(*,*)
        write(*,*) ' Max abs error = ', emax
        write(*,*) ' R.m.s.  error = ', erms
        write(*,*)
    endif


    call print_result_dp('Test inversion (lapinv0)', emax, atol=7.0e-7)

    deallocate(vorx)
    deallocate(potx)
    deallocate(potxe)
    deallocate(ds)

end program test_lapinv0_2
