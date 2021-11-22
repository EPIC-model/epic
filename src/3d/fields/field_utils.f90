module field_utils
    use parameters, only : dx, dxi, extent, lower, nx, ny, nz
    use inversion_utils, only : fftxyp2s, fftxys2p, diffx, diffy
    use constants, only : f12
    implicit none

    contains

        subroutine gradient(f, fgrad)
            double precision, intent(in)  :: f(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: fgrad(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision              :: fs(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1)

            ! copy field (we use the output temporarily)
            fgrad(:, :, :, 1) = f(0:nz, :, :)

            ! compute spectral field:
            call fftxyp2s(fgrad(:, :, :, 1), fs)

            ! df/dx
            call diffx(fs, ds)                      ! df/dx in spectral space
            call fftxys2p(ds, fgrad(:, :, :, 1))    ! df/dx in physical space

            ! df/dy
            call diffy(fs, ds)                      ! df/dy in spectral space
            call fftxys2p(ds, fgrad(:, :, :, 2))    ! df/dy in physical space

            ! df/dz (central difference)
            fgrad(:, :, :, 3) = f12 * dxi(3) * (f(1:nz+1, :, :) - f(-1:nz-1, :, :))
        end subroutine gradient

        ! Get the lower index of the cell the parcel is in.
        ! This subroutine does not take x periodicity into account.
        ! @param[in] pos position of the parcel
        ! @param[out] i lower, zonal cell index
        ! @param[out] j lower, vertical cell index
        subroutine get_index(pos, i, j, k)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: i, j, k
            integer                       :: idx(3)

            idx = floor((pos - lower) * dxi)

            i = idx(1)
            j = idx(2)
            k = idx(3)
        end subroutine get_index


        ! Do periodic shift of the index
        ! @param[inout] ii zonal grid point indices
        ! @param[inout] jj meridional grid point indices
        subroutine periodic_index_shift(ii, jj)
            integer, intent(inout) :: ii(:), jj(:)

            ! account for x / y periodicity:
            ! -1          --> nx-1 / ny-1
            !  0          --> 0
            ! nx+1 / ny+1 --> 1
            ! nx / ny     --> 0
            ! nx-1 / ny-1 --> nx-1 / ny-1
            ii = mod(ii + nx, nx)
            jj = mod(jj + ny, ny)

        end subroutine periodic_index_shift


        ! Get the coordinate of a grid point (i, j, k).
        ! @param[in] i zonal cell index
        ! @param[in] j meridional cell index
        ! @param[in] k vertical cell index
        ! @param[out] pos position of (i, j, k) in the domain
        subroutine get_position(i, j, k, pos)
            integer,          intent(in)  :: i, j, k
            double precision, intent(out) :: pos(3)

            pos = lower + dble((/i, j, k/)) * dx

        end subroutine get_position

end module field_utils
