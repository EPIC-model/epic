module inversion_mod
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use inversion_utils
    use parameters, only : nx, ny, nz, dxi, l_bndry_zeta_zero
    use physics, only : f_cor
    use constants, only : zero, two, f12
    use sta2dfft, only : dct, dst
    use sta3dfft, only : rkz, rkzi, ztrig, zfactors, diffx, diffy, fftxyp2s, fftxys2p
    use fields
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given the vorticity vector field (svor) in spectral space, this
        ! returns the associated velocity field (velog) as well as vorticity
        ! and the velocity gradient tensor (velgradg) in physical space (vortg)
        ! Note: the vorticity is modified to be solenoidal and spectrally filtered.
        subroutine vor2vel
            double precision :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))          ! semi-spectral
            double precision :: bs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))          ! semi-spectral
            double precision :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))          ! semi-spectral
            double precision :: es(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))          ! semi-spectral
            double precision :: cs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))          ! semi-spectral
            double precision :: svel(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), n_dim) ! semi-spectral
            double precision :: svor(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), n_dim) ! mixed spectral
            double precision :: ubar(0:nz), vbar(0:nz)
            integer          :: iz, nc, kx, ky, kz

            call start_timer(vor2vel_timer)

            !----------------------------------------------------------
            ! Decompose initial vorticity and filter spectrally:
            do nc = 1, n_dim
                call field_decompose_physical(vortg(:, :, :, nc), svor(:, :, :, nc))
                svor(:, :, :, nc) = filt * svor(:, :, :, nc)
            enddo

            !----------------------------------------------------------
            ! Enforce solenoidality
            ! A, B, C are vorticities
            ! D = B_x - A_y; E = C_z
            ! A = k2l2i * (E_x + D_y) and B = k2l2i * (E_y - D_x) --> A_x + B_y + C_z = zero
            call diffx(svor(:, :, :, I_Y), as) ! as = B_x
            call diffy(svor(:, :, :, I_X), bs) ! bs = A_y
            !$omp parallel workshare
            ds = as - bs                     ! ds = D
            cs = svor(:, :, :, I_Z)
            !$omp end parallel workshare
            call field_combine_semi_spectral(cs)
            call central_diffz(cs, es)                     ! es = E
            call field_decompose_semi_spectral(es)

            ! ubar and vbar are used here to store the mean x and y components of the vorticity
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                ubar = svor(:, 0, 0, I_X)
                vbar = svor(:, 0, 0, I_Y)
            endif

            call diffx(es, svor(:, :, :, I_X)) ! E_x
            call diffy(ds, cs)               ! cs = D_y
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
               svor(iz, :, :, I_X) = k2l2i * (svor(iz, :, :, I_X) + cs(iz, :, :))
            enddo
            !$omp end parallel do

            call diffy(es, svor(:, :, :, I_Y)) ! E_y
            call diffx(ds, cs)                 ! D_x

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
               svor(iz, :, :, I_Y) = k2l2i * (svor(iz, :, :, I_Y) - cs(iz, :, :))
            enddo
            !$omp end parallel do

            ! bring back the mean x and y components of the vorticity
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                svor(:, 0, 0, I_X) = ubar
                svor(:, 0, 0, I_Y) = vbar
            endif

            !----------------------------------------------------------
            ! Combine vorticity in physical space:
            do nc = 1, n_dim
                call field_combine_physical(svor(:, :, :, nc), vortg(:, :, :, nc))
            enddo

            !----------------------------------------------------------
            !Form source term for inversion of vertical velocity -> ds:
            call diffy(svor(:, :, :, I_X), ds)
            call diffx(svor(:, :, :, I_Y), es)
            !$omp parallel workshare
            ds = ds - es
            !$omp end parallel workshare

            !Calculate the boundary contributions of the source to the vertical velocity (bs)
            !and its derivative (es) in semi-spectral space:
            !$omp parallel do private(iz)  default(shared)
            do iz = 1, nz-1
                bs(iz, :, :) = ds(0, :, :) *  thetam(iz, :, :) + ds(nz, :, :) *  thetap(iz, :, :)
            enddo
            !$omp end parallel do

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                es(iz, :, :) = ds(0, :, :) * dthetam(iz, :, :) + ds(nz, :, :) * dthetap(iz, :, :)
            enddo
            !$omp end parallel do

            !Invert Laplacian to find the part of w expressible as a sine series:
            !$omp parallel workshare
            ds(1:nz-1, :, :) = green(1:nz-1, :, :) * ds(1:nz-1, :, :)
            !$omp end parallel workshare

            ! Calculate d/dz of this sine series:
            !$omp parallel workshare
            as(0, :, :) = zero
            !$omp end parallel workshare
            !$omp parallel do private(iz)  default(shared)
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * ds(kz, :, :)
            enddo
            !$omp end parallel do
            !$omp parallel workshare
            as(nz, :, :) = zero
            !$omp end parallel workshare

            !FFT these quantities back to semi-spectral space:
            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, as(0:nz, ky, kx), ztrig, zfactors)
                    call dst(1, nz, ds(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            ! Combine vertical velocity (ds) and its derivative (es) given the sine and linear parts:
            !$omp parallel workshare
            ds(0     , :, :) = zero
            ds(1:nz-1, :, :) = ds(1:nz-1, :, :) + bs(1:nz-1, :, :)
            ds(nz    , :, :) = zero
            es = es + as

            ! Get complete zeta field in semi-spectral space
            cs = svor(:, :, :, I_Z)
            !$omp end parallel workshare
            call field_combine_semi_spectral(cs)

            !----------------------------------------------------------------------
            !Define horizontally-averaged flow by integrating the horizontal vorticity:

            !First integrate the sine series in svor(1:nz-1, 0, 0, I_X & I_Y):
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                ubar(0) = zero
                vbar(0) = zero
                ubar(1:nz-1) = -rkzi * svor(1:nz-1, 0, 0, I_Y)
                vbar(1:nz-1) =  rkzi * svor(1:nz-1, 0, 0, I_X)
                ubar(nz) = zero
                vbar(nz) = zero

                !Transform to semi-spectral space as a cosine series:
                call dct(1, nz, ubar, ztrig, zfactors)
                call dct(1, nz, vbar, ztrig, zfactors)

                !Add contribution from the linear function connecting the boundary values:
                ubar = ubar + svor(nz, 0, 0, I_Y) * gamtop - svor(0, 0, 0, I_Y) * gambot
                vbar = vbar - svor(nz, 0, 0, I_X) * gamtop + svor(0, 0, 0, I_X) * gambot
            endif

            !-------------------------------------------------------
            !Find x velocity component "u":
            call diffx(es, as)
            call diffy(cs, bs)

            !$omp parallel do private(iz) default(shared)
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) + bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                as(:, 0, 0) = ubar
            endif

            !Store spectral form of "u":
            !$omp parallel workshare
            svel(:, :, :, I_X) = as
            !$omp end parallel workshare

            !Get "u" in physical space:
            call fftxys2p(as, velog(:, :, :, I_X))

            !-------------------------------------------------------
            !Find y velocity component "v":
            call diffy(es, as)
            call diffx(cs, bs)

            !$omp parallel do private(iz) default(shared)
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) - bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                as(:, 0, 0) = vbar
            endif

            !Store spectral form of "v":
            !$omp parallel workshare
            svel(:, :, :, I_Y) = as
            !$omp end parallel workshare

            !Get "v" in physical space:
            call fftxys2p(as, velog(:, :, :, I_Y))

            !-------------------------------------------------------
            !Store spectral form of "w":
            !$omp parallel workshare
            svel(:, :, :, I_Z) = ds
            !$omp end parallel workshare

            !Get "w" in physical space:
            call fftxys2p(ds, velog(:, :, :, I_Z))

            !=================================================================================

            ! compute the velocity gradient tensor
            call vel2vgrad(svel)

            ! use extrapolation in u and v and anti-symmetry in w to fill z grid points outside domain:
            !$omp parallel workshare
            velog(-1, :, :, I_X) =  two * velog(0, :, :, I_X) - velog(1, :, :, I_X) ! u
            velog(-1, :, :, I_Y) =  two * velog(0, :, :, I_Y) - velog(1, :, :, I_Y) ! v
            velog(-1, :, :, I_Z) = -velog(1, :, :, I_Z) ! w
            velog(nz+1, :, :, I_X) = two * velog(nz, :, :, I_X) - velog(nz-1, :, :, I_X) ! u
            velog(nz+1, :, :, I_Y) = two * velog(nz, :, :, I_Y) - velog(nz-1, :, :, I_Y) ! v
            velog(nz+1, :, :, I_Z) = -velog(nz-1, :, :, I_Z) ! w
            !$omp end parallel workshare

            call stop_timer(vor2vel_timer)

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded velocity gradient tensor
        subroutine vel2vgrad(svel)
            double precision, intent(in) :: svel(box%lo(3):box%hi(3), & ! 0:nz
                                                 box%lo(2):box%hi(2), &
                                                 box%lo(1):box%hi(1), &
                                                 n_dim) ! velocity in semi-spectral space
            double precision             :: ds(box%lo(3):box%hi(3), & ! 0:nz
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1)) ! semi-spectral derivatives

            ! x component:
            call diffx(svel(:, :, :, I_X), ds)           ! u_x = du/dx in semi-spectral space
            call fftxys2p(ds, velgradg(:, :, :, I_DUDX)) ! u_x in physical space

            call diffy(svel(:, :, :, I_X), ds)              ! u_y = du/dy in semi-spectral space
            call fftxys2p(ds, velgradg(:, :, :, I_DUDY)) ! u_y in physical space

            call diffx(svel(:, :, :, I_Z), ds)              ! w_x = dw/dx in semi-spectral space
            call fftxys2p(ds, velgradg(:, :, :, I_DWDX)) ! w_x in physical space

            ! use extrapolation in du/dx and du/dy to fill z grid points outside domain:
            !$omp parallel workshare
            velgradg(  -1, :, :, I_DUDX) =  two * velgradg( 0, :, :, I_DUDX) - velgradg(   1, :, :, I_DUDX)
            velgradg(nz+1, :, :, I_DUDX) =  two * velgradg(nz, :, :, I_DUDX) - velgradg(nz-1, :, :, I_DUDX)
            velgradg(  -1, :, :, I_DUDY) =  two * velgradg( 0, :, :, I_DUDY) - velgradg(   1, :, :, I_DUDY)
            velgradg(nz+1, :, :, I_DUDY) =  two * velgradg(nz, :, :, I_DUDY) - velgradg(nz-1, :, :, I_DUDY)

            ! use anti-symmetry for dw/dx to fill z grid points outside domain:
            velgradg(  -1, :, :, I_DWDX) = -velgradg(   1, :, :, I_DWDX)
            velgradg(nz+1, :, :, I_DWDX) = -velgradg(nz-1, :, :, I_DWDX)
            !$omp end parallel workshare

            ! y & z components:
            call diffy(svel(:, :, :, I_Y), ds)              ! v_y = dv/dy in semi-spectral space
            call fftxys2p(ds, velgradg(:, :, :, I_DVDY)) ! v_y in physical space

            call diffy(svel(:, :, :, I_Z), ds)              ! w_y = dw/dy in semi-spectral space
            call fftxys2p(ds, velgradg(:, :, :, I_DWDY)) ! w_y in physical space

            !$omp parallel workshare
            ! use extrapolation in dv/dy to fill z grid points outside domain:
            velgradg(  -1, :, :, I_DVDY) = two * velgradg( 0, :, :, I_DVDY) - velgradg(   1, :, :, I_DVDY)
            velgradg(nz+1, :, :, I_DVDY) = two * velgradg(nz, :, :, I_DVDY) - velgradg(nz-1, :, :, I_DVDY)

            ! use anti-symmetry in dw/dy to fill z grid points outside domain:
            ! w_y(-1) = -w_y(1) and w_y(nz+1) = -w_y(nz-1)
            velgradg(  -1, :, :, I_DWDY) = -velgradg(   1, :, :, I_DWDY)
            velgradg(nz+1, :, :, I_DWDY) = -velgradg(nz-1, :, :, I_DWDY)
            !$omp end parallel workshare

        end subroutine vel2vgrad

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine vorticity_tendency
            double precision :: f(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), n_dim)

            call start_timer(vtend_timer)

            ! Eqs. 10 and 11 of MPIC paper
            f(:, : , :, I_X) = (vortg(:, :, :, I_X) + f_cor(I_X)) * velog(:, :, :, I_X)
            f(:, : , :, I_Y) = (vortg(:, :, :, I_Y) + f_cor(I_Y)) * velog(:, :, :, I_X) + tbuoyg
            f(:, : , :, I_Z) = (vortg(:, :, :, I_Z) + f_cor(I_Z)) * velog(:, :, :, I_X)

            call divergence(f, vtend(0:nz, :, :, I_X))

            f(:, : , :, I_X) = (vortg(:, :, :, I_X) + f_cor(I_X)) * velog(:, :, :, I_Y) - tbuoyg
            f(:, : , :, I_Y) = (vortg(:, :, :, I_Y) + f_cor(I_Y)) * velog(:, :, :, I_Y)
            f(:, : , :, I_Z) = (vortg(:, :, :, I_Z) + f_cor(I_Z)) * velog(:, :, :, I_Y)

           call divergence(f, vtend(0:nz, :, :, I_Y))

            f(:, : , :, I_X) = (vortg(:, :, :, I_X) + f_cor(I_X)) * velog(:, :, :, I_Z)
            f(:, : , :, I_Y) = (vortg(:, :, :, I_Y) + f_cor(I_Y)) * velog(:, :, :, I_Z)
            f(:, : , :, I_Z) = (vortg(:, :, :, I_Z) + f_cor(I_Z)) * velog(:, :, :, I_Z)

            call divergence(f, vtend(0:nz, :, :, I_Z))

            !-------------------------------------------------------
            ! Set dzeta/dt = 0 on the boundary if required:
            if (l_bndry_zeta_zero(1)) then
                vtend(0, :, :, I_Z) = zero
            endif

            if (l_bndry_zeta_zero(2)) then
                vtend(nz, :, :, I_Z) = zero
            endif

            !-------------------------------------------------------
            ! Extrapolate to halo grid points:
            !$omp parallel workshare
            vtend(-1,   :, :, :) = two * vtend(0,  :, :, :) - vtend(1,    :, :, :)
            vtend(nz+1, :, :, :) = two * vtend(nz, :, :, :) - vtend(nz-1, :, :, :)
            !$omp end parallel workshare
        end subroutine vorticity_tendency

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Note: f is overwritten
        subroutine divergence(f, div)
            double precision, intent(inout) :: f(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), n_dim)
            double precision, intent(out)   :: div(0:nz, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision                :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            ! calculate df1/dx
            call fftxyp2s(f(:, :, :, I_X), fs)
            call diffx(fs, ds)
            call fftxys2p(ds, f(:, :, :, I_X))

            ! calculate df2/dy
            call fftxyp2s(f(:, :, :, I_Y), fs)
            call diffy(fs, ds)
            call fftxys2p(ds, f(:, :, :, I_Y))

            ! calculate df3/dz
            call central_diffz(f(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), I_Z), &
                             div(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            ! div = df1/dx + df2/dy + df3/dz
            div = f(0:nz, :, :, I_X) + f(0:nz, :, :, I_Y) + div

          end subroutine divergence

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a divergent flow field (ud, vd, wd) = grad(phi) where
        ! Lap(phi) = div (given).
        subroutine diverge(div,  ud, vd, wd)
            double precision, intent(inout)  :: div(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision, intent(out)    :: ud(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)), &
                                                vd(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)), &
                                                wd(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision                 :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                 :: us(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                                                vs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                                                ws(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            !------------------------------------------------------------------
            ! Convert phi to spectral space (in x & y) as ds:
            call fftxyp2s(div, ds)

            if ((1 >= box%lo(1)) .and. (1 <= box%hi(1)) .and. &
                (1 >= box%lo(2)) .and. (1 <= box%hi(2))) then
                ds(:, 1, 1) = zero
            endif

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

            ! Compute z derivative by central differences:
            call central_diffz(ds, ws)

            ! Set vertical boundary values to zero
            ws(0,  :, :) = zero
            ws(nz, :, :) = zero

            ! Reverse FFT to define z velocity component wd:
            call fftxys2p(ws, wd)

        end subroutine

end module inversion_mod
