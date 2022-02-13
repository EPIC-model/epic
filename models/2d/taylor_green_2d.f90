! =============================================================================
!                           2D TAYLOR-GREEN FLOW
!
!                   u(x, y) = A * cos(ax + d) * sin(by + e)
!                   v(x, y) = B * sin(ax + d) * cos(by + e)
!
!                   The code uses following variable names:
!                       tg_flow%freq  = (/a, b/)
!                       tg_flow%amp   = (/A, B/)
!                       tg_flow%phase = (/d, e/)
! =============================================================================
module taylor_green_2d
    use netcdf_writer
    use constants, only : pi, f12, zero, one, two
    implicit none

    private
        type flow_type
            double precision :: amp(2) = (/f12, -one/)      ! amplitudes
            double precision :: freq(2) = (/two, one/)      ! frequencies
            double precision :: phase(2) = (/zero, zero/)   ! phase shift
        end type flow_type

        type(flow_type) :: tg_flow

        integer :: vor_id

        double precision, parameter :: hpi = f12 * pi

    public :: get_flow_velocity,    &
              get_flow_gradient,    &
              get_flow_vorticity,   &
              taylor_green_init,    &
              tg_flow


    contains
        subroutine taylor_green_init(ncfname, ncid, dimids, nx, nz, origin, dx)
            character(*),     intent(in)    :: ncfname
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, nz
            double precision, intent(in)    :: origin(2)
            double precision, intent(in)    :: dx(2)
            double precision                :: pos(2)
            double precision                :: vortg(0:nz, 0:nx-1)
            integer                         :: i, j

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='vorticity',                    &
                                       long_name='',                        &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=vor_id)

            do j = 0, nz
                do i = 0, nx - 1
                    pos = origin + dx * dble((/i, j/))
                    vortg(j, i) = get_flow_vorticity(pos)
                enddo
            enddo

            call write_netcdf_dataset(ncid, vor_id, vortg)

        end subroutine taylor_green_init


        function get_flow_velocity(pos) result(vel)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: vel(2)

            call get_flow_pos(pos, xx, zz)

            vel(1) = tg_flow%amp(1) * dcos(xx) * dsin(zz)
            vel(2) = tg_flow%amp(2) * dsin(xx) * dcos(zz)
        end function get_flow_velocity

        ! grad ordering : dudx, dudy, dvdx, dvdy
        function get_flow_gradient(pos) result(grad)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: grad(4)

            call get_flow_pos(pos, xx, zz)

            ! du/dx = - a * A * sin(xx) * sin(zz)
            grad(1) = - tg_flow%freq(1) * tg_flow%amp(1) * dsin(xx) * dsin(zz)

            ! du/dy = b * A * cos(xx) * cos(zz)
            grad(2) = tg_flow%freq(2) * tg_flow%amp(1) * dcos(xx) * dcos(zz)

            ! dv/dx = a * B * cos(xx) * np.cos(zz)
            grad(3) = tg_flow%freq(1) * tg_flow%amp(2) * dcos(xx) * dcos(zz)

            ! dv/dy = - b * B * sin(xx) * sin(zz)
            grad(4) = - tg_flow%freq(2) * tg_flow%amp(2) * dsin(xx) * dsin(zz)

        end function get_flow_gradient

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: omega

            call get_flow_pos(pos, xx, zz)

            omega = (tg_flow%amp(2) * tg_flow%freq(1)     &
                   - tg_flow%amp(1) * tg_flow%freq(2))    &
                   * dcos(xx) * dcos(zz)
        end function get_flow_vorticity

        subroutine get_flow_pos(pos, xx, zz)
            double precision, intent(in) :: pos(2)
            double precision, intent(out) :: xx, zz

            xx = tg_flow%freq(1) * pos(1) + tg_flow%phase(1) + hpi
            zz = tg_flow%freq(2) * pos(2) + tg_flow%phase(2)

        end subroutine get_flow_pos

end module taylor_green_2d
