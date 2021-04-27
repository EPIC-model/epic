! =============================================================================
!             This module contains various interpolation methods.
! =============================================================================
module interpl_methods
    use parameters, only : dx, dxi
    use fields, only : get_index                &
                     , get_position             &
                     , periodic_index_shift
    implicit none

    contains

        !
        ! tri-linear interpolation
        !
        subroutine trilinear(pos, ji, weight, ngp)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ji(2, 4)
            double precision, intent(out) :: weight(4)
            integer,          intent(out) :: ngp
            double precision              :: xy(2)

            ! (j, i)
            ji(:, 1) = get_index(pos)
            xy = get_position(ji(:, 1))
            weight(1) = product(1.0 - abs(pos - xy) * dxi)

            ! (j, i+1)
            ji(:, 2) = ji(:, 1)
            ji(2, 2) = ji(2, 1) + 1
            xy = get_position(ji(:, 2))
            weight(2) = product(1.0 - abs(pos - xy) * dxi)

            ! (j+1, i)
            ji(:, 3) = ji(:, 1)
            ji(1, 3) = ji(1, 3) + 1
            xy = get_position(ji(:, 3))
            weight(3) = product(1.0 - abs(pos - xy) * dxi)

            ! (j+1, i+1)
            ji(:, 4) = ji(:, 2)
            ji(1, 4) = ji(1, 4) + 1
            xy = get_position(ji(:, 4))
            weight(4) = product(1.0 - abs(pos - xy) * dxi)

            ngp = 4

            ! account for x periodicity
            call periodic_index_shift(ji, ngp)

        end subroutine trilinear

end module interpl_methods
