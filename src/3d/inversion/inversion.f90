module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz
    use phys_parameters, only : ft_cor, f_cor
    use constants, only : zero, two, three, four
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given the vorticity vector field (vortg) in physical space, this
        ! returns the associated velocity field (velog) and the velocity
        ! gradient tensor (velgradg).  Note: the
        ! vorticity is modified to be solenoidal and spectrally filtered.
        subroutine vor2vel(vortg,  velog,  velgradg)
            double precision, intent(inout)  :: vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3)
!             ox(0:nz, ny, nx), oy(0:nz, ny, nx), oz(0:nz, ny, nx)
            double precision, intent(out) :: velog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
!             uu(0:nz, ny, nx), vv(0:nz, ny, nx), ww(0:nz, ny, nx)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
!             us(0:nz, nx, ny), vs(0:nz, nx, ny), ws(0:nz, nx, ny)
            double precision              :: svelog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision              :: as(0:nz, nx, ny), bs(0:nz, nx, ny), cs(0:nz, nx, ny)
            double precision              :: ds(0:nz, nx, ny), es(0:nz, nx, ny), fs(0:nz, nx, ny)
            double precision              :: astop(nx, ny), bstop(nx, ny)
            double precision              :: asbot(nx, ny), bsbot(nx, ny)
            double precision              :: ubar(0:nz), vbar(0:nz)
            integer                       :: iz

            call start_timer(vor2vel_timer)

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

            ! compute the velocity gradient tensor
            call vel2vgrad(vortg, svelog, velgradg)

            call stop_timer(vor2vel_timer)

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded velocity gradient tensor
        subroutine vel2vgrad(vortg, svelog, velgradg)
            double precision, intent(in)  :: vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in)  :: svelog(-1:nz+1, 0:nx-1, 0:ny-1, 3)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives

            ! x component:
            call diffx(svelog(0:nz, :, :, 1), ds)      ! u_x = du/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 1)) ! u_x in physical space

            call diffy(svelog(0:nz, :, :, 1), ds)      ! u_y = du/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 2)) ! u_y in physical space

            call diffx(svelog(0:nz, :, :, 3), ds)      ! w_x = dw/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 4)) ! w_x in physical space

            ! use symmetry to fill z grid points outside domain:
            ! u_x(-1) =  u_x(1) and u_x(nz+1) =  u_x(nz-1)
            ! u_y(-1) =  u_y(1) and u_y(nz+1) =  u_y(nz-1)
            ! w_x(-1) = -w_x(1) and w_x(nz+1) = -w_x(nz-1)
            velgradg(  -1, :, :, 1) = velgradg(   1, :, :, 1) ! lower boundary du/dx
            velgradg(nz+1, :, :, 1) = velgradg(nz-1, :, :, 1) ! upper boundary du/dx
            velgradg(  -1, :, :, 2) = velgradg(   1, :, :, 2) ! lower boundary du/dy
            velgradg(nz+1, :, :, 2) = velgradg(nz-1, :, :, 2) ! upper boundary du/dy
            velgradg(  -1, :, :, 4) = velgradg(   1, :, :, 4) ! lower boundary dw/dx
            velgradg(nz+1, :, :, 4) = velgradg(nz-1, :, :, 4) ! upper boundary dw/dx


            ! y & z components:

            call diffy(svelog(0:nz, :, :, 2), ds)      ! v_y = dv/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 3)) ! v_y in physical space

            call diffy(svelog(0:nz, :, :, 3), ds)      ! w_y = dw/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 5)) ! w_y in physical space

            ! use symmetry to fill z grid points outside domain:
            ! v_y(-1) =  v_y(1) and v_y(nz+1) =  v_y(nz-1)
            ! w_y(-1) = -w_y(1) and w_y(nz+1) = -w_y(nz-1)
            velgradg(  -1, :, :, 3) = velgradg(   1, :, :, 3) ! lower boundary dv/dy
            velgradg(nz+1, :, :, 3) = velgradg(nz-1, :, :, 3) ! upper boundary dv/dy
            velgradg(  -1, :, :, 5) = velgradg(   1, :, :, 5) ! lower boundary dw/dy
            velgradg(nz+1, :, :, 5) = velgradg(nz-1, :, :, 5) ! upper boundary dw/dy

        end subroutine vel2vgrad

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency:
        subroutine vorticity_tendency(vortg, tbuoyg, velgradg, vtend)
            double precision, intent(in)  :: vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in)  :: tbuoyg(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision              :: b(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision, intent(out) :: vtend(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision              :: bs(0:nz, 0:nx-1, 0:ny-1) ! spectral buoyancy
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives
            double precision              :: db(0:nz, 0:ny-1, 0:nx-1) ! buoyancy derivatives

            call start_timer(vtend_timer)

            ! copy buoyancy
            b = tbuoyg(0:nz, :, :)

            ! Compute spectral buoyancy (bs):
            call fftxyp2s(b, bs)

            call diffy(bs, ds)                      ! b_y = db/dy in spectral space
            call fftxys2p(ds, db)                   ! db = b_y in physical space

            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, 1) =  vortg(0:nz, :, :, 1)           * velgradg(0:nz, :, :, 1) & ! \omegax * du/dx
                                 + (vortg(0:nz, :, :, 2) + ft_cor) * velgradg(0:nz, :, :, 2) & ! \omegay * du/dy
                                 + (vortg(0:nz, :, :, 3) +  f_cor) *                         &
                                            (vortg(0:nz, :, :, 2) + velgradg(0:nz, :, :, 4)) & ! \omegaz * du/dz
                                 + db                                                          ! db/dy
            !$omp end workshare
            !$omp end parallel

            call diffx(bs, ds)                      ! b_x = db/dx in spectral space
            call fftxys2p(ds, db)                   ! db = b_x in physical space

            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, 2) =  vortg(0:nz, :, :, 1)           *                         &
                                            (vortg(0:nz, :, :, 3) + velgradg(0:nz, :, :, 2)) & ! \omegax * dv/dx
                                 + (vortg(0:nz, :, :, 2) + ft_cor) * velgradg(0:nz, :, :, 3) & ! \omegay * dv/dy
                                 + (vortg(0:nz, :, :, 3) + f_cor)  *                         &
                                            (velgradg(0:nz, :, :, 5) - vortg(0:nz, :, :, 1)) & ! \omegaz * dv/dz
                                 - db                                                          ! dbdx

            vtend(0:nz, :, :, 3) =  vortg(0:nz, :, :, 1)           * velgradg(0:nz, :, :, 4) & ! \omegax * dw/dx
                                 + (vortg(0:nz, :, :, 2) + ft_cor) * velgradg(0:nz, :, :, 5) & ! \omegay * dw/dy
                                 - (vortg(0:nz, :, :, 3) + f_cor)  *                         &
                                         (velgradg(0:nz, :, :, 1) + velgradg(0:nz, :, :, 3))   ! \omegaz * dw/dz

            !$omp end workshare
            !$omp end parallel

            ! Extrapolate to halo grid points
            vtend(-1,   :, :, :) = two * vtend(0,  :, :, :) - vtend(1,    :, :, :)
            vtend(nz+1, :, :, :) = two * vtend(nz, :, :, :) - vtend(nz-1, :, :, :)

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

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
