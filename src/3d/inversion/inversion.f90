module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz
    use constants, only : zero, two, f12
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,   &
               db_timer

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given the vorticity vector field (vortg) in physical space, this
        ! returns the associated velocity field (velog) and the velocity
        ! gradient tensor (velgradg).
        subroutine vor2vel(vortg,  velog,  velgradg)
            double precision, intent(in)    :: vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out)   :: velog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out)   :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision                :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: as(0:nz, 0:nx-1, 0:ny-1) &
                                             , bs(0:nz, 0:nx-1, 0:ny-1) &
                                             , cs(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ds(0:nz, 0:nx-1, 0:ny-1) &
                                             , es(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: dstop(0:nx-1, 0:ny-1), dsbot(0:nx-1, 0:ny-1)
            double precision                :: ubar(0:nz), vbar(0:nz)
            double precision                :: uavg, vavg
            integer                         :: iz

            call start_timer(vor2vel_timer)

            ! Copy vorticity to velocity field to perform FFT transforms
            ! (FFT transforms overwrite the input array)
            velog = vortg

            !------------------------------------------------------------------
            !Convert vorticity to spectral space as (as, bs, cs): (velog is overwritten in this operation)
            call fftxyp2s(velog(0:nz, :, :, 1), as)
            call fftxyp2s(velog(0:nz, :, :, 2), bs)
            call fftxyp2s(velog(0:nz, :, :, 3), cs)

            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (bs(iz, 0, 0) + bs(iz+1, 0, 0))
                vbar(iz+1) = vbar(iz) - dz2 * (as(iz, 0, 0) + as(iz+1, 0, 0))
            enddo

            ! remove the mean value to have zero net momentum
            uavg = sum(ubar(1:nz-1) + f12 * ubar(nz)) / dble(nz)
            vavg = sum(vbar(1:nz-1) + f12 * vbar(nz)) / dble(nz)
            do iz = 0, nz
                ubar(iz) = ubar(iz) - uavg
                vbar(iz) = vbar(iz) - vavg
            enddo

            !Form source term for inversion of vertical velocity:
            call diffy(as, ds)
            call diffx(bs, es)

            !$omp parallel
            !$omp workshare
            ds = ds - es
            !$omp end workshare
            !$omp end parallel

            !as & bs are now free to re-use

            !Save boundary values for z derivative of w below:
            dsbot = ds(0,  :, :)
            dstop = ds(nz, :, :)

            !Invert to find vertical velocity \hat{w} (store in ds, spectrally):
            call lapinv0(ds)

            !Find \hat{w}' (store in es, spectrally):
            call diffz0(ds, es, dsbot, dstop)

            !Find x velocity component \hat{u}:
            call diffx(es, as)
            call diffy(cs, bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) + bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = ubar

            svelog(:, :, :, 1) = as

            !Get u in physical space:
            call fftxys2p(as, velog(0:nz, :, :, 1))

            !Find y velocity component \hat{v}:
            call diffy(es, as)
            call diffx(cs, bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) - bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = vbar

            svelog(:, :, :, 2) = as

            !Get v in physical space:
            call fftxys2p(as, velog(0:nz, :, :, 2))

            svelog(:, :, :, 3) = ds

            !Get w in physical space:
            call fftxys2p(ds, velog(0:nz, :, :, 3))

            ! compute the velocity gradient tensor
            call vel2vgrad(svelog, velgradg)

            ! use extrapolation in u and v and anti-symmetry in w to fill z grid points outside domain:
            velog(-1, :, :, 1) =  two * velog(0, :, :, 1) - velog(1, :, :, 1) ! u
            velog(-1, :, :, 2) =  two * velog(0, :, :, 2) - velog(1, :, :, 2) ! v
            velog(-1, :, :, 3) = -velog(1, :, :, 3) ! w
            velog(nz+1, :, :, 1) = two * velog(nz, :, :, 1) - velog(nz-1, :, :, 1) ! u
            velog(nz+1, :, :, 2) = two * velog(nz, :, :, 2) - velog(nz-1, :, :, 2) ! v
            velog(nz+1, :, :, 3) = -velog(nz-1, :, :, 3) ! w

            call stop_timer(vor2vel_timer)

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded velocity gradient tensor
        subroutine vel2vgrad(svelog, velgradg)
            double precision, intent(in)  :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives

            ! x component:
            call diffx(svelog(:, :, :, 1), ds)         ! u_x = du/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 1)) ! u_x in physical space

            call diffy(svelog(:, :, :, 1), ds)         ! u_y = du/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 2)) ! u_y in physical space

            call diffx(svelog(:, :, :, 3), ds)         ! w_x = dw/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 4)) ! w_x in physical space

            ! use extrapolation in du/dx and du/dy to fill z grid points outside domain:
            velgradg(  -1, :, :, 1) =  two * velgradg( 0, :, :, 1) - velgradg(   1, :, :, 1) ! lower boundary du/dx
            velgradg(nz+1, :, :, 1) =  two * velgradg(nz, :, :, 1) - velgradg(nz-1, :, :, 1) ! upper boundary du/dx
            velgradg(  -1, :, :, 2) =  two * velgradg( 0, :, :, 2) - velgradg(   1, :, :, 2) ! lower boundary du/dy
            velgradg(nz+1, :, :, 2) =  two * velgradg(nz, :, :, 2) - velgradg(nz-1, :, :, 2) ! upper boundary du/dy

            ! use anti-symmetry for dw/dx to fill z grid points outside domain:
            velgradg(  -1, :, :, 4) = -velgradg(   1, :, :, 4) ! lower boundary dw/dx
            velgradg(nz+1, :, :, 4) = -velgradg(nz-1, :, :, 4) ! upper boundary dw/dx

            ! y & z components:
            call diffy(svelog(:, :, :, 2), ds)         ! v_y = dv/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 3)) ! v_y in physical space

            call diffy(svelog(:, :, :, 3), ds)         ! w_y = dw/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 5)) ! w_y in physical space

            ! use extrapolation in dv/dy to fill z grid points outside domain:
            velgradg(  -1, :, :, 3) = two * velgradg( 0, :, :, 3) - velgradg(   1, :, :, 3) ! lower boundary dv/dy
            velgradg(nz+1, :, :, 3) = two * velgradg(nz, :, :, 3) - velgradg(nz-1, :, :, 3) ! upper boundary dv/dy

            ! use anti-symmetry in dw/dy to fill z grid points outside domain:
            ! w_y(-1) = -w_y(1) and w_y(nz+1) = -w_y(nz-1)
            velgradg(  -1, :, :, 5) = -velgradg(   1, :, :, 5) ! lower boundary dw/dy
            velgradg(nz+1, :, :, 5) = -velgradg(nz-1, :, :, 5) ! upper boundary dw/dy

        end subroutine vel2vgrad

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded buoyancy derivatives db/dx and db/dy
        subroutine buoyancy_derivatives(tbuoyg, dbdx, dbdy)
            double precision, intent(in)  :: tbuoyg(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: dbdx(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: dbdy(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision              :: b(0:nz, 0:ny-1, 0:nx-1)    ! buoyancy in physical space
            double precision              :: bs(0:nz, 0:nx-1, 0:ny-1)   ! buoyancy in spectral space
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1)   ! buoyancy derivative in spectral space

            call start_timer(db_timer)

            ! copy buoyancy
            b = tbuoyg(0:nz, :, :)

            ! Compute spectral buoyancy (bs):
            call fftxyp2s(b, bs)

            call diffy(bs, ds)                      ! b_y = db/dy in spectral space
            call fftxys2p(ds, dbdy(0:nz, :, :))     ! db = b_y in physical space

            call diffx(bs, ds)                      ! b_x = db/dx in spectral space
            call fftxys2p(ds, dbdx(0:nz, :, :))     ! db = b_x in physical space

            ! Extrapolate to halo grid points
            dbdy(-1,   :, :) = two * dbdy(0,  :, :) - dbdy(1,    :, :)
            dbdy(nz+1, :, :) = two * dbdy(nz, :, :) - dbdy(nz-1, :, :)

            dbdx(-1,   :, :) = two * dbdx(0,  :, :) - dbdx(1,    :, :)
            dbdx(nz+1, :, :) = two * dbdx(nz, :, :) - dbdx(nz-1, :, :)

            call stop_timer(db_timer)

        end subroutine buoyancy_derivatives

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
