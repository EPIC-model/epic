! =============================================================================
!                       Straka ("Cold bubble") test case
!
!   This module sets up a cold bubble as described in Straka et al (1993).
!
!   Straka, J.M., Wilhelmson, R.B., Wicker, L.J., Anderson, J.R. and
!   Droegemeier, K.K. (1993), Numerical solutions of a non-linear density
!   current: A benchmark solution and comparisons. Int. J. Numer. Meth. Fluids,
!   17: 1-22. https://doi.org/10.1002/fld.1650170103
!
! =============================================================================

module straka
    use phys_constants
    use constants
    use h5_writer
    implicit none

    private

    type flow_type
        double precision :: theta_ref = 300.0d0               ![Kelvin] reference potential temperature
        double precision :: theta_max = 15.0d0                ![Kelvin] max. pot. temp. perturbation
        double precision :: center(2) = (/zero, 3000.0d0/)    ![m] sphere center (x, z)
        double precision :: radii(2)  = (/4000.0d0, 2000.d0/) ![m] ellipse radii (x, z)
    end type flow_type

    type(flow_type) :: straka_flow

    public :: straka_init, straka_flow

    contains

        subroutine straka_init(h5fname, nx, nz, origin, dx)
            character(*),     intent(in) :: h5fname
            integer,          intent(in) :: nx, nz
            double precision, intent(in) :: origin(2)
            double precision, intent(in) :: dx(2)
            double precision             :: pos(2)
            double precision             :: xc, xr, zc, zr, L
            double precision             :: dtheta, theta_ref, theta_max
            double precision             :: buoyg(0:nz, 0:nx-1)
            integer                      :: i, j
            integer(hid_t)               :: h5handle
            logical                      :: exists = .true.

            ! check whether file exists
            inquire(file=h5fname, exist=exists)
            if (exists) then
                print *, "File '" // h5fname // "'already exists."
                stop
            endif

            ! in metres
            xc = straka_flow%center(1)
            xr = straka_flow%radii(2)
            zc = straka_flow%center(2)
            zr = straka_flow%radii(2)

            ![Kelvin] reference potential temperature
            theta_ref = straka_flow%theta_ref

            theta_max = -0.5d0 * straka_flow%theta_max

            do j = 0, nz
                do i = 0, nx - 1
                    pos = origin + dx * dble((/i, j/))

                    L = ((pos(1) - xc) / xr) ** 2 &
                      + ((pos(2) - zc) / zr) ** 2

                    L = dsqrt(L)

                    ! potential temperature perturbation
                    dtheta = zero

                    if (L <= one) then
                        dtheta = theta_max * (dcos(pi * L) + one)
                    endif

                    ! MPIC paper:
                    ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                    ! (dtheta = theta - theta_ref)
                    buoyg(j, i) = g_c * dtheta / theta_ref
                enddo
            enddo

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5handle)
            call write_h5_dataset_2d(h5handle, '/', 'buoyancy', buoyg)
            call close_h5_file(h5handle)

        end subroutine straka_init
end module
