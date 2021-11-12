module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz
    use constants, only : zero, three, four
    implicit none

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given the vorticity vector field (vortg) in physical space, this
        ! returns the associated velocity field (velog) and a spectral
        ! copy (svelog),  which is additionally filtered.  Note: the
        ! vorticity is modified to be solenoidal and spectrally filtered.
        subroutine vor2vel(vortg,  velog,  svelog)
            double precision, intent(inout)  :: vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3)
!             ox(0:nz, ny, nx), oy(0:nz, ny, nx), oz(0:nz, ny, nx)
            double precision, intent(out) :: velog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
!             uu(0:nz, ny, nx), vv(0:nz, ny, nx), ww(0:nz, ny, nx)
            double precision, intent(out) :: svelog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
!             us(0:nz, nx, ny), vs(0:nz, nx, ny), ws(0:nz, nx, ny)
            double precision              :: as(0:nz, nx, ny), bs(0:nz, nx, ny), cs(0:nz, nx, ny)
            double precision              :: ds(0:nz, nx, ny), es(0:nz, nx, ny), fs(0:nz, nx, ny)
            double precision              :: astop(nx, ny), bstop(nx, ny)
            double precision              :: asbot(nx, ny), bsbot(nx, ny)
            double precision              :: ubar(0:nz), vbar(0:nz)
            integer                       :: iz

            !------------------------------------------------------------------
            !Convert vorticity to spectral space as (as, bs, cs):
            call fftxyp2s(vortg(0:nz, :, :, 1), as)
            call fftxyp2s(vortg(0:nz, :, :, 2), bs)
            call fftxyp2s(vortg(0:nz, :, :, 3), cs)

            !Add -grad(lambda) where Laplace(lambda) = div(vortg) to
            !enforce the solenoidal condition on the vorticity field:
            call diffx(as, ds)
            call diffy(bs, es)

            !For the vertical parcel vorticity, use 2nd-order accurate
            !differencing with extrapolation at the boundaries:
            !$omp parallel shared(fs, cs) private(iz)
            !$omp do
            do iz = 1, nz-1
                fs(iz, :, :) = hdzi * (cs(iz+1, :, :) - cs(iz-1, :, :))
            enddo
            !$omp end do
            !$omp end parallel

            fs(0, :, :) = hdzi * (four * cs(1, :, :) - cs(2, :, :) - three * cs(0, :, :))
            fs(nz, :, :) = hdzi * (cs(nz-2, :, :) + three * cs(nz, :, :) - four * cs(nz-1, :, :))

            !Form div(vortg):
            !$omp parallel
            !$omp workshare
            fs = ds + es + fs
            !$omp end workshare
            !$omp end parallel

            !Remove horizontally-averaged part (plays no role):
            fs(:, 1, 1) = zero

            !Invert Lap(lambda) = div(vortg) assuming dlambda/dz = 0 at the
            !boundaries (store solution lambda in fs):
            call lapinv1(fs)

            !Filter lambda:
            !$omp parallel shared(fs, filt, nz) private(iz)  default(none)
            !$omp do
            do iz = 0, nz
                fs(iz, :, :) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Subtract grad(lambda) to enforce div(vortg) = 0:
            call diffx(fs, ds)
            !$omp parallel
            !$omp workshare
            as = as - ds
            !$omp end workshare
            !$omp end parallel

            call diffy(fs, ds)
            !$omp parallel
            !$omp workshare
            bs = bs - ds
            !$omp end workshare
            !$omp end parallel

            call diffz1(fs, ds)
            !$omp parallel
            !$omp workshare
            cs = cs - ds
            !$omp end workshare
            !$omp end parallel
            !Ensure horizontal average of vertical vorticity is zero:
            cs(:, 1, 1) = zero

            !Compute spectrally filtered vorticity in physical space:
            !$omp parallel shared(ds, es, fs, as, bs, cs, filt, nz) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                ds(iz, :, :) = filt * as(iz, :, :)
                es(iz, :, :) = filt * bs(iz, :, :)
                fs(iz, :, :) = filt * cs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Save boundary values of x and y vorticity for z derivatives of A & B:
            asbot = as(0, :, :)
            bsbot = bs(0, :, :)
            astop = as(nz, :, :)
            bstop = bs(nz, :, :)

            !Return corrected vorticity to physical space:
            call fftxys2p(ds, vortg(0:nz, :, :, 1))
            call fftxys2p(es, vortg(0:nz, :, :, 2))
            call fftxys2p(fs, vortg(0:nz, :, :, 3))

            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (es(iz, 1, 1) + es(iz+1, 1, 1))
                vbar(iz+1) = vbar(iz) - dz2 * (ds(iz, 1, 1) + ds(iz+1, 1, 1))
            enddo

            !-----------------------------------------------------------------
            !Invert vorticity to find vector potential (A, B, C) -> (as, bs, cs):
            call lapinv0(as)
            call lapinv0(bs)
            call lapinv1(cs)

            !------------------------------------------------------------
            !Compute x velocity component, u = B_z - C_y:
            call diffy(cs, ds)
            call diffz0(bs, es, bsbot, bstop)
            !bsbot & bstop contain spectral y vorticity component at z_min and z_max
            !$omp parallel
            !$omp workshare
            fs = es - ds
            !$omp end workshare
            !$omp end parallel
            !Add horizontally-averaged flow:
            fs(:, 1, 1) = ubar
            !$omp parallel shared(svelog, fs, filt, nz) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                svelog(iz, :, :, 1) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get u in physical space:
            call fftxys2p(fs, velog(0:nz, :, :, 1))

            !------------------------------------------------------------
            !Compute y velocity component, v = C_x - A_z:
            call diffx(cs, ds)
            call diffz0(as, es, asbot, astop)
            !asbot & astop contain spectral x vorticity component at z_min and z_max
            !$omp parallel
            !$omp workshare
            fs = ds - es
            !$omp end workshare
            !$omp end parallel
            !Add horizontally-averaged flow:
            fs(:, 1, 1) = vbar
            !$omp parallel shared(svelog, fs, filt, nz) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                svelog(iz, :, :, 2) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get v in physical space:
            call fftxys2p(fs, velog(0:nz, :, :, 2))

            !------------------------------------------------------------
            !Compute z velocity component, w = A_y - B_x:
            call diffx(bs, ds)
            call diffy(as, es)
            !$omp parallel
            !$omp workshare
            fs = es - ds
            !$omp end workshare
            !$omp end parallel

            !$omp parallel shared(svelog, fs, filt, nz) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                svelog(iz, :, :, 3) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get w in physical space:
            call fftxys2p(fs, velog(0:nz, :, :, 3))
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a divergent flow field (ud, vd, wd) = grad(phi) where
        ! Lap(phi) = div (given).
        subroutine diverge(div,  ud, vd, wd)
            double precision, intent(inout)  :: div(0:nz, ny, nx)
            double precision, intent(out)    :: ud(0:nz, ny, nx), vd(0:nz, ny, nx), wd(0:nz, ny, nx)
            double precision                 :: ds(0:nz, nx, ny)
            double precision                 :: us(0:nz, nx, ny), vs(0:nz, nx, ny), ws(0:nz, nx, ny)
            double precision                 :: wbar(0:nz)

            !------------------------------------------------------------------
            ! Convert phi to spectral space (in x & y) as ds:
            call fftxyp2s(div, ds)

            ! Compute the x & y-independent part of ds by integration:
            call vertint(ds(:, 1, 1), wbar)

            ! Invert Laplace's operator semi-spectrally with compact differences:
            call lapinv1(ds)

            ! Compute x derivative spectrally:
            call diffx(ds, us)

            ! Reverse FFT to define x velocity component ud:
            call fftxys2p(us, ud)

            ! Compute y derivative spectrally:
            call diffy(ds, vs)

            ! Reverse FFT to define y velocity component vd:
            call fftxys2p(vs, vd)

            ! Compute z derivative by compact differences:
            call diffz1(ds, ws)

            ! Add on the x and y-independent part of wd:
            ws(:, 1, 1) = ws(:, 1, 1) + wbar

            ! Reverse FFT to define z velocity component wd:
            call fftxys2p(ws, wd)

        end subroutine

end module inversion_mod
