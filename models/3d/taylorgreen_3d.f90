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
module taylorgreen_3d
    use h5_writer
    use constants, only : pi, f12, zero, one, two
    implicit none

    private
        type flow_type
            double precision :: amp(3) = (/one, one, -two/)      ! amplitudes
            double precision :: freq(3) = (/one, one, one/)      ! frequencies
            double precision :: phase(3) = (/zero, zero, zero/)  ! phase shift
        end type flow_type

        type(flow_type) :: tg_flow


        double precision, parameter :: hpi = f12 * pi

    public :: get_flow_vorticity,   &
              taylorgreen_init,     &
              tg_flow


    contains
        subroutine taylorgreen_init(h5handle, nx, ny, nz, origin, dx)
            integer(hid_t),   intent(inout) :: h5handle
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            double precision                :: pos(3)
            double precision                :: vortg(0:nz, 0:ny-1, 0:nx-1, 3)
            integer                         :: i, j, k

            do i = 0, nx - 1
                do j = 0, ny - 1
                    do k = 0, nz
                        pos = origin + dx * dble((/i, j, k/))
                        vortg(k, j, i, :) = get_flow_vorticity(pos)
                    enddo
                enddo
            enddo

            call write_h5_dataset(h5handle, '/', 'vorticity', vortg)

        end subroutine taylorgreen_init

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

end module taylorgreen_3d
