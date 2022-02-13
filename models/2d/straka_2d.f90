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

module straka_2d
    use phys_constants
    use constants
    use netcdf_writer
    implicit none

    private

    type flow_type
        double precision :: dtheta_max = 15.0d0               ![K] max. pot. temp. perturbation
        double precision :: center(2) = (/zero, 3000.0d0/)    ![m] sphere center (x, z)
        double precision :: radii(2)  = (/4000.0d0, 2000.d0/) ![m] ellipse radii (x, z)
    end type flow_type

    integer :: buo_id

    type(flow_type) :: straka_flow

    public :: straka_init, straka_flow

    contains

        subroutine straka_init(ncid, dimids, nx, nz, origin, dx)
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, nz
            double precision, intent(in)    :: origin(2)
            double precision, intent(in)    :: dx(2)
            double precision                :: pos(2)
            double precision                :: xc, xr, zc, zr, L
            double precision                :: dtheta, dtheta_max
            double precision                :: buoyg(0:nz, 0:nx-1)
            integer                         :: i, j

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='buoyancy',                &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=buo_id)

            ! in metres
            xc = straka_flow%center(1)
            xr = straka_flow%radii(1)
            zc = straka_flow%center(2)
            zr = straka_flow%radii(2)

            dtheta_max = -f12 * straka_flow%dtheta_max

            do j = 0, nz
                do i = 0, nx - 1
                    pos = origin + dx * dble((/i, j/))

                    L = ((pos(1) - xc) / xr) ** 2 &
                      + ((pos(2) - zc) / zr) ** 2

                    L = dsqrt(L)

                    ! potential temperature perturbation
                    dtheta = zero

                    if (L <= one) then
                        dtheta = dtheta_max * (dcos(pi * L) + one)
                    endif

                    ! MPIC paper:
                    ! liquid-water buoyancy is defined by b = g * (theta − theta_l0) / theta_l0
                    ! (dtheta = theta - theta_l0)
                    buoyg(j, i) = gravity * dtheta / theta_l0
                enddo
            enddo

            call write_netcdf_dataset(ncid, buo_id, buoyg)

        end subroutine straka_init
end module
