! ====================================================================================
!                           3D TAYLOR-GREEN FLOW
!
! Taylor-Green flow:
!     u(x, y, z) = A * cos(a * x + d) * sin(b * y + e) * sin(c * z + f)
!     v(x, y, z) = B * sin(a * x + d) * cos(b * y + e) * sin(c * z + f)
!     w(x, y, z) = C * sin(a * x + d) * sin(b * y + e) * cos(c * z + f)
! The vorticity of this flow is
!    xi(x, y, z) = (b * C - c * B) * sin(a * x + d) * cos(b * y + e) * cos(c * z + f)
!   eta(x, y, z) = (c * A - a * C) * cos(a * x + d) * sin(b * y + e) * cos(c * z + f)
!  zeta(x, y, z) = (a * B - b * A) * cos(a * x + d) * cos(b * y + e) * sin(c * z + f)
!
!                   The code uses following variable names:
!                       tg_flow%freq  = (/a, b, c/)
!                       tg_flow%amp   = (/A, B, C/)
!                       tg_flow%phase = (/d, e, f/)
! ====================================================================================
module taylor_green_3d
    use netcdf_writer
    use constants, only : pi, f12, zero, one, two
    implicit none

    private
        type flow_type
            double precision :: amp(3) = (/one, one, -two/)      ! amplitudes
            double precision :: freq(3) = (/one, one, one/)      ! frequencies
            double precision :: phase(3) = (/zero, zero, zero/)  ! phase shift
        end type flow_type

        type(flow_type) :: tg_flow

        integer :: x_vor_id, y_vor_id, z_vor_id


        double precision, parameter :: hpi = f12 * pi

    public :: get_flow_vorticity,   &
              taylor_green_init,    &
              tg_flow


    contains
        subroutine taylor_green_init(ncfname, ncid, dimids, nx, ny, nz, origin, dx)
            character(*),     intent(in)    :: ncfname
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            double precision                :: pos(3)
            double precision                :: vortg(0:nz, 0:ny-1, 0:nx-1, 3)
            integer                         :: i, j, k

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_vorticity',                  &
                                       long_name='x vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_vorticity',                  &
                                       long_name='y vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_vorticity',                  &
                                       long_name='z vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vor_id)

            call close_definition(ncid)

            do i = 0, nx - 1
                do j = 0, ny - 1
                    do k = 0, nz
                        pos = origin + dx * dble((/i, j, k/))
                        vortg(k, j, i, :) = get_flow_vorticity(pos)
                    enddo
                enddo
            enddo

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            call write_netcdf_dataset(ncid, x_vor_id, vortg(:, :, :, 1))
            call write_netcdf_dataset(ncid, y_vor_id, vortg(:, :, :, 2))
            call write_netcdf_dataset(ncid, z_vor_id, vortg(:, :, :, 3))

        end subroutine taylor_green_init

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(3)
            double precision             :: xx, yy, zz
            double precision             :: omega(3)
            double precision             :: A, B, C

            call get_flow_pos(pos, xx, yy, zz)

            A = tg_flow%freq(2) * tg_flow%amp(3) &
              - tg_flow%freq(3) * tg_flow%amp(2)

            B = tg_flow%freq(3) * tg_flow%amp(1) &
              - tg_flow%freq(1) * tg_flow%amp(3)

            C = tg_flow%freq(1) * tg_flow%amp(2) &
              - tg_flow%freq(2) * tg_flow%amp(1)

            omega(1) = A * dsin(xx) * dcos(yy) * dcos(zz)
            omega(2) = B * dcos(xx) * dsin(yy) * dcos(zz)
            omega(3) = C * dcos(xx) * dcos(yy) * dsin(zz)

        end function get_flow_vorticity

        subroutine get_flow_pos(pos, xx, yy, zz)
            double precision, intent(in) :: pos(3)
            double precision, intent(out) :: xx, yy, zz

            xx = tg_flow%freq(1) * pos(1) + tg_flow%phase(1) + hpi
            yy = tg_flow%freq(2) * pos(2) + tg_flow%phase(2) + hpi
            zz = tg_flow%freq(3) * pos(3) + tg_flow%phase(3)

        end subroutine get_flow_pos

end module taylor_green_3d
