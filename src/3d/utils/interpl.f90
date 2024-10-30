!==============================================================================
!                   Module for interpolation methods
!==============================================================================
module interpl
    use constants, only : one
    use parameters, only : lower, dxi
    implicit none

    contains

        ! Tri-linear interpolation
        ! @param[in] pos position vector
        ! @param[out] ii zonal lower grid point for interoplation
        ! @param[out] jj meridional lower grid point for interpolation
        ! @param[out] kk vertical lower grid point for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine trilinear(pos, ii, jj, kk, ww)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: ii, jj, kk
            double precision, intent(out) :: ww(0:1, 0:1, 0:1)
            double precision              :: xyz(3)
            double precision              :: w00, w10, w01, w11
            double precision              :: px, py, pz, pxc, pyc, pzc


            ! (i, j, k)
            xyz = (pos - lower) * dxi
            ii = floor(xyz(1))
            jj = floor(xyz(2))
            kk = floor(xyz(3))

            px = xyz(1) - dble(ii)
            pxc = one - px

            py = xyz(2) - dble(jj)
            pyc = one - py

            pz = xyz(3) - dble(kk)
            pzc = one - pz

            w00 = pyc * pxc
            w10 = pyc * px
            w01 = py * pxc
            w11 = py * px

            ! Note order of indices is k,j,i
            ww(0, 0, 0) = pzc * w00
            ww(0, 0, 1) = pzc * w10
            ww(0, 1, 0) = pzc * w01
            ww(0, 1, 1) = pzc * w11
            ww(1, 0, 0) = pz * w00
            ww(1, 0, 1) = pz * w10
            ww(1, 1, 0) = pz * w01
            ww(1, 1, 1) = pz * w11

        end subroutine trilinear

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Bi-linear interpolation
        ! @param[in] pos position vector
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj meridional grid points for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine bilinear(pos, ii, jj, ww)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ii, jj
            double precision, intent(out) :: ww(0:1, 0:1)
            double precision              :: xy(2)
            double precision              :: px, py, pxc, pyc


            ! (i, j)
            xy = (pos - lower(1:2)) * dxi(1:2)
            ii = floor(xy(1))
            jj = floor(xy(2))

            px = xy(1) - dble(ii)
            pxc = one - px

            py = xy(2) - dble(jj)
            pyc = one - py

            ! Note order of indices is j,i
            ww(0, 0) = pyc * pxc
            ww(0, 1) = pyc * px
            ww(1, 0) = py  * pxc
            ww(1, 1) = py  * px

        end subroutine bilinear

end module interpl
