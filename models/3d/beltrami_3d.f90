! ====================================================================================
!                           3D BELTRAMI FLOW
!
! Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) + l*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
! The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
!
!                   The code uses following variable names:
!                       beltrami_flow%alpha = alpha
!                       beltrami_flow%k     = k
!                       beltrami_flow%l     = l
!                       beltrami_flow%m     = m
! ====================================================================================
module beltrami_3d
    use netcdf_writer
    use constants, only : pi, f12, zero, one, two
    implicit none

    private
        type beltrami_type
            integer          :: k, l, m
            double precision :: alpha
        end type beltrami_type

        type(beltrami_type) :: beltrami_flow

        double precision :: fk2l2, kk, ll, mm

        integer :: x_vor_id, y_vor_id, z_vor_id

        double precision, parameter :: hpi = f12 * pi

    public :: get_flow_vorticity,   &
              beltrami_init,        &
              beltrami_flow


    contains
        subroutine beltrami_init(ncid, dimids, nx, ny, nz, origin, dx)
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


            kk = dble(beltrami_flow%k)
            ll = dble(beltrami_flow%l)
            mm = dble(beltrami_flow%m)
            fk2l2 = one / dble(beltrami_flow%k ** 2 + beltrami_flow%l ** 2)

            do i = 0, nx - 1
                do j = 0, ny - 1
                    do k = 0, nz
                        pos = origin + dx * dble((/i, j, k/))
                        vortg(k, j, i, :) = get_flow_vorticity(pos)
                    enddo
                enddo
            enddo

            call write_netcdf_dataset(ncid, x_vor_id, vortg(:, :, :, 1))
            call write_netcdf_dataset(ncid, y_vor_id, vortg(:, :, :, 2))
            call write_netcdf_dataset(ncid, z_vor_id, vortg(:, :, :, 3))

        end subroutine beltrami_init

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(3)
            double precision             :: xx, yy, zz
            double precision             :: omega(3)
            double precision             :: cosmz, sinmz, sinkxly, coskxly

            call get_flow_pos(pos, xx, yy, zz)

            cosmz = dcos(mm * zz)
            sinmz = dsin(mm * zz)
            sinkxly = dsin(kk * xx + ll * yy)
            coskxly = dcos(kk * xx + ll * yy)

            omega(1) = fk2l2 * (kk * mm * sinmz - ll * beltrami_flow%alpha * cosmz) * sinkxly
            omega(2) = fk2l2 * (kk * mm * sinmz + ll * beltrami_flow%alpha * cosmz) * sinkxly
            omega(3) = cosmz * coskxly

        end function get_flow_vorticity

        subroutine get_flow_pos(pos, xx, yy, zz)
            double precision, intent(in) :: pos(3)
            double precision, intent(out) :: xx, yy, zz

            xx = pos(1) + hpi
            yy = pos(2) + hpi
            zz = pos(3) + hpi

        end subroutine get_flow_pos

end module beltrami_3d
