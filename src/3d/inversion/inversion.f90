module inversion_mod
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use inversion_utils
    use parameters, only : nx, ny, nz, dxi
    use physics, only : f_cor
    use constants, only : zero, two, f12
    use sta2dfft, only : dct, dst
    use fields
    use options, only : l_flux
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
            double precision :: as(0:nz, 0:nx-1, 0:ny-1)        ! semi-spectral
            double precision :: bs(0:nz, 0:nx-1, 0:ny-1)        ! semi-spectral
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)        ! semi-spectral
            double precision :: es(0:nz, 0:nx-1, 0:ny-1)        ! semi-spectral
            double precision :: cs(0:nz, 0:nx-1, 0:ny-1)        ! semi-spectral
            double precision :: ubar(0:nz), vbar(0:nz)
            double precision :: svel(0:nz, 0:nx-1, 0:ny-1, n_dim) ! velocity in semi-spectral space
            double precision :: svor(0:nz, 0:nx-1, 0:ny-1, n_dim) ! vorticity in mixed spectral space
            integer          :: iz, nc, kx, ky, kz

            call start_timer(vor2vel_timer)

            !----------------------------------------------------------
            ! Decompose initial vorticity and filter spectrally:
            do nc = 1, n_dim
                call field_decompose_physical(vortg(0:nz, :, :, nc), svor(:, :, :, nc))
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
            !call field_combine_semi_spectral(cs)
            !call diffz(cs, es)                     ! es = E
            call spectral_diffz(cs, es)
            !call field_decompose_semi_spectral(es)

            ! ubar and vbar are used here to store the mean x and y components of the vorticity
            ubar = svor(:, 0, 0, I_X)
            vbar = svor(:, 0, 0, I_Y)

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
            svor(:, 0, 0, I_X) = ubar
            svor(:, 0, 0, I_Y) = vbar

            !----------------------------------------------------------
            ! Combine vorticity in physical space:
            do nc = 1, n_dim
                call field_combine_physical(svor(:, :, :, nc), vortg(0:nz, :, :, nc))
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
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dct(1, nz, as(0:nz, kx, ky), ztrig, zfactors)
                    call dst(1, nz, ds(1:nz, kx, ky), ztrig, zfactors)
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
            as(:, 0, 0) = ubar

            !Store spectral form of "u":
            !$omp parallel workshare
            svel(:, :, :, I_X) = as
            !$omp end parallel workshare

            !Get "u" in physical space:
            call fftxys2p(as, velog(0:nz, :, :, I_X))

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
            as(:, 0, 0) = vbar

            !Store spectral form of "v":
            !$omp parallel workshare
            svel(:, :, :, I_Y) = as
            !$omp end parallel workshare

            !Get "v" in physical space:
            call fftxys2p(as, velog(0:nz, :, :, I_Y))

            !-------------------------------------------------------
            !Store spectral form of "w":
            !$omp parallel workshare
            svel(:, :, :, I_Z) = ds
            !$omp end parallel workshare

            !Get "w" in physical space:
            call fftxys2p(ds, velog(0:nz, :, :, I_Z))

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
            double precision, intent(in) :: svel(0:nz, 0:nx-1, 0:ny-1, n_dim) ! velocity in semi-spectral space
            double precision             :: ds(0:nz, 0:nx-1, 0:ny-1)          ! semi-spectral derivatives

            ! x component:
            call diffx(svel(:, :, :, I_X), ds)           ! u_x = du/dx in semi-spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, I_DUDX)) ! u_x in physical space

            call diffy(svel(:, :, :, I_X), ds)              ! u_y = du/dy in semi-spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, I_DUDY)) ! u_y in physical space

            call diffx(svel(:, :, :, I_Z), ds)              ! w_x = dw/dx in semi-spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, I_DWDX)) ! w_x in physical space

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
            call fftxys2p(ds, velgradg(0:nz, :, :, I_DVDY)) ! v_y in physical space

            call diffy(svel(:, :, :, I_Z), ds)              ! w_y = dw/dy in semi-spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, I_DWDY)) ! w_y in physical space

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

          if (l_flux) then
             call vorticity_tendency_flux
          else
             call vorticity_tendency_ogradu_diffz
!             call vorticity_tendency_ogradu
          endif
        end subroutine vorticity_tendency


        subroutine vorticity_tendency_ogradu_diffz
            double precision :: b(0:nz, 0:ny-1, 0:nx-1)
            double precision :: bs(0:nz, 0:nx-1, 0:ny-1) ! spectral buoyancy
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives
            double precision :: db(0:nz, 0:ny-1, 0:nx-1) ! buoyancy derivatives
            double precision :: dds(0:nz, 0:nx-1, 0:ny-1)
            double precision :: ddf(0:nz, 0:ny-1, 0:nx-1)
            call start_timer(vtend_timer)

            ! copy buoyancy
            b = tbuoyg(0:nz, :, :)
            
            ! Compute spectral buoyancy (bs):
            call fftxyp2s(b, bs)
            call diffy(bs, ds)                      ! b_y = db/dy in spectral space
            call fftxys2p(ds, db)                   ! db = b_y in physical space

            ! du/dz
            call field_decompose_physical(velog(0:nz, :, :, 1), ds)
            call spectral_diffz(ds, dds)
            call field_combine_physical(dds, ddf) ! ddf = du/dz

!            ! du/dz = \omegay + dw/dx
!            ddf = vortg(0:nz, :, :, 2) + velgradg(0:nz, :, :, 4)

            !$omp parallel 
            !$omp workshare
            vtend(0:nz, :, :, 1) =  vortg(0:nz, :, :, 1)             * velgradg(0:nz, :, :, 1) & ! \omegax * du/dx
                                 + (vortg(0:nz, :, :, 2) + f_cor(2)) * velgradg(0:nz, :, :, 2) & ! \omegay * du/dy
                                 + (vortg(0:nz, :, :, 3) + f_cor(3)) * ddf           !          & ! \omegaz * du/dz
!                                 + db                                                            ! db/dy
            !$omp end workshare
            !$omp end parallel
            call diffx(bs, ds)                      ! b_x = db/dx in spectral space
            call fftxys2p(ds, db)                   ! db = b_x in physical space


            ! dv/dz
            call field_decompose_physical(velog(0:nz, :, :, 2), ds)
            call spectral_diffz(ds, dds)
            call field_combine_physical(dds, ddf) ! ddf = dv/dz

            ! dv/dz = dw/dy - \omegax
 !           ddf = velgradg(0:nz, :, :, 5) - vortg(0:nz, :, :, 1)
            
            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, 2) =  vortg(0:nz, :, :, 1)             * (vortg(0:nz, :, :, 3) + velgradg(0:nz, :, :, 2)) & ! \omegax * dv/dx (dv/dx = du/dy + \omegaz)
                                 + (vortg(0:nz, :, :, 2) + f_cor(2)) * velgradg(0:nz, :, :, 3) & ! \omegay * dv/dy
                                 + (vortg(0:nz, :, :, 3) + f_cor(3)) * ddf                   !  & ! \omegaz * dv/dz
                                ! - db
            !$omp end workshare
            !$omp end parallel

            ! dw/dz
            call field_decompose_physical(velog(0:nz, :, :, 3), ds)
            call spectral_diffz(ds, dds)
            call field_combine_physical(dds, ddf) ! ddf = dw/dz

            ! dw/dz = - du/dx - dv/dy
!            ddf = - velgradg(0:nz, :, :, 1) - velgradg(0:nz, :, :, 3)

            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, 3) =  vortg(0:nz, :, :, 1)             * velgradg(0:nz, :, :, 4) & ! \omegax * dw/dx
                                 + (vortg(0:nz, :, :, 2) + f_cor(2)) * velgradg(0:nz, :, :, 5) & ! \omegay * dw/dy
                                 + (vortg(0:nz, :, :, 3) + f_cor(3)) * ddf                       ! \omegaz * dw/dz (dw/dz = - du/dx - dv/dy)
            !$omp end workshare
            !$omp end parallel
            
            ! Extrapolate to halo grid points
            vtend(-1,   :, :, :) = two * vtend(0,  :, :, :) - vtend(1,    :, :, :)
            vtend(nz+1, :, :, :) = two * vtend(nz, :, :, :) - vtend(nz-1, :, :, :)
            
            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency_ogradu_diffz

        subroutine vorticity_tendency_ogradu
            double precision :: b(0:nz, 0:ny-1, 0:nx-1)
            double precision :: bs(0:nz, 0:nx-1, 0:ny-1) ! spectral buoyancy
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives
            double precision :: db(0:nz, 0:ny-1, 0:nx-1) ! buoyancy derivatives

            call start_timer(vtend_timer)

            ! copy buoyancy
            b = tbuoyg(0:nz, :, :)

            ! Compute spectral buoyancy (bs):
            call fftxyp2s(b, bs)

            call diffy(bs, ds)    ! b_y = db/dy in spectral space
            call fftxys2p(ds, db) ! db = b_y in physical space

            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, I_X) =  vortg(0:nz, :, :, I_X) * velgradg(0:nz, :, :, I_DUDX)    & ! \xi * du/dx
                                   + (vortg(0:nz, :, :, I_Y) + f_cor(I_Y))                     &
                                        * (vortg(0:nz, :, :, I_Z) + velgradg(0:nz, :, :, I_Y)) & ! \eta * dv/dx
                                   + (vortg(0:nz, :, :, I_Z) + f_cor(I_Z))                     &
                                        * velgradg(0:nz, :, :, I_DWDX)                         & ! \zeta * dw/dx
                                   + db                                                          ! db/dy
            !$omp end workshare
            !$omp end parallel

            call diffx(bs, ds)    ! b_x = db/dx in spectral space
            call fftxys2p(ds, db) ! db = b_x in physical space

            !$omp parallel
            !$omp workshare
            vtend(0:nz, :, :, I_Y) =  vortg(0:nz, :, :, I_X) * velgradg(0:nz, :, :, I_Y) & ! \xi * du/dy
                                   + (vortg(0:nz, :, :, I_Y) + f_cor(I_Y))               &
                                        * velgradg(0:nz, :, :, I_DVDY)                   & ! \eta * dv/dy
                                   + (vortg(0:nz, :, :, I_Z) + f_cor(I_Z))               &
                                        * velgradg(0:nz, :, :, I_DWDY)                   & ! \zeta * dw/dy
                                   - db                                                    ! dbdx

            vtend(0:nz, :, :, I_Z) =  vortg(0:nz, :, :, I_X) * velgradg(0:nz, :, :, I_DWDX) & ! \xi * dw/dx
                                   + (vortg(0:nz, :, :, I_Y) + f_cor(I_Y))                  &
                                        * velgradg(0:nz, :, :, I_DWDY)                      & ! \eta * dw/dy
                                   - (vortg(0:nz, :, :, I_Z) + f_cor(I_Z))  *               &
                                    (velgradg(0:nz, :, :, I_DUDX) + velgradg(0:nz, :, :, I_DVDY)) ! \zeta * dw/dz

            !$omp end workshare
            !$omp end parallel

            ! Extrapolate to halo grid points
            vtend(-1,   :, :, :) = two * vtend(0,  :, :, :) - vtend(1,    :, :, :)
            vtend(nz+1, :, :, :) = two * vtend(nz, :, :, :) - vtend(nz-1, :, :, :)

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency_ogradu

        subroutine vorticity_tendency_flux
            double precision :: f(-1:nz+1, 0:ny-1, 0:nx-1, 3)

            call start_timer(vtend_timer)

            ! Eqs. 10 and 11 of MPIC paper
            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 1)
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 1) + tbuoyg
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 1)

            call divergence(f, vtend(0:nz, :, :, 1))

            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 2) - tbuoyg
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 2)
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 2)

           call divergence(f, vtend(0:nz, :, :, 2))

            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 3)
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 3)
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 3)

            call divergence(f, vtend(0:nz, :, :, 3))

            !-------------------------------------------------------
            ! Extrapolate to halo grid points:
            !$omp parallel workshare
            vtend(-1,   :, :, :) = two * vtend(0,  :, :, :) - vtend(1,    :, :, :)
            vtend(nz+1, :, :, :) = two * vtend(nz, :, :, :) - vtend(nz-1, :, :, :)
            !$omp end parallel workshare
        end subroutine vorticity_tendency_flux

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine divergence(f, div)
            double precision, intent(in)  :: f(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out) :: div(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: df(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: dds(0:nz, 0:nx-1, 0:ny-1)
!            integer                       :: i

!            ! calculate df/dx with central differencing
!            do i = 1, nx-2
!                div(0:nz, 0:ny-1, i) = f12 * dxi(1) * (f(0:nz, 0:ny-1, i+1, 1) - f(0:nz, 0:ny-1, i-1, 1))
!            enddo
!            div(0:nz, 0:ny-1, 0)    = f12 * dxi(1) * (f(0:nz, 0:ny-1, 1, 1) - f(0:nz, 0:ny-1, nx-1, 1))
!            div(0:nz, 0:ny-1, nx-1) = f12 * dxi(1) * (f(0:nz, 0:ny-1, 0, 1) - f(0:nz, 0:ny-1, nx-2, 1))


            df = f(0:nz, :, :, 1)
            call fftxyp2s(df, ds)
            call diffx(ds, dds)
            call fftxys2p(dds, div)

            df = f(0:nz, :, :, 2)
            call fftxyp2s(df, ds)
            call diffy(ds, dds)
            call fftxys2p(dds, df)
            div = div + df
           
!            ! calculate df/dy with central differencing
!!            do i = 1, ny-2
!                df(0:nz, i, 0:nx-1) = f12 * dxi(2) * (f(0:nz, i+1, 0:nx-1, 2) - f(0:nz, i-1, 0:nx-1, 2))
!            enddo
!            df(0:nz, 0,    0:nx-1) = f12 * dxi(2) * (f(0:nz, 1, 0:nx-1, 2) - f(0:nz, ny-1, 0:nx-1, 2))
!            df(0:nz, ny-1, 0:nx-1) = f12 * dxi(2) * (f(0:nz, 0, 0:nx-1, 2) - f(0:nz, ny-2, 0:nx-1, 2))
!            div = div + df

            ! calculate df/dz with central differencing
!!            call diffz(f(0:nz, :, :, 3), df)
!            do i = 0, nz
!                df(i, 0:ny-1, 0:nx-1) = f12 * dxi(3) * (f(i+1, 0:ny-1, 0:nx-1, 3) - f(i-1, 0:ny-1, 0:nx-1, 3))
            !            enddo

            call field_decompose_physical(f(0:nz, :, :, 3), ds)
            call spectral_diffz(ds, dds)
            call field_combine_physical(dds, df)

            div = div + df

          end subroutine divergence


          subroutine spectral_diffz(fs, ds)
!            double precision, intent(in)  :: f(0:nz, 0:ny-1, 0:nx-1)
!            double precision, intent(out) :: df(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(in)  :: fs(0:nz, 0:nx-1, 0:ny-1) ! f in mixed-spectral space
            double precision, intent(out) :: ds(0:nz, 0:nx-1, 0:ny-1) ! derivative linear part
            double precision              :: as(0:nz, 0:nx-1, 0:ny-1) ! derivative sine-part
            integer :: kx, ky, kz, iz

!            call field_decompose_physical(f, fs)

            !Calculate the boundary contributions of the derivative (ds) in semi-spectral space:
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                ds(iz, :, :) = fs(0, :, :) * dthetam(iz, :, :) + fs(nz, :, :) * dthetap(iz, :, :)
            enddo
            !$omp end parallel do
            
            ! Calculate d/dz of this sine series:
            !$omp parallel workshare
            as(0, :, :) = zero
            !$omp end parallel workshare
            !$omp parallel do private(kz)  default(shared)
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * fs(kz, :, :)
            enddo
            !$omp end parallel do
            !$omp parallel workshare
            as(nz, :, :) = zero
            !$omp end parallel workshare

            !FFT these quantities back to semi-spectral space:
            !$omp parallel do collapse(2) private(kx, ky)
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dct(1, nz, as(0:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do
            
            ! Combine vertical derivative (es) given the sine and linear parts:
            !omp parallel workshare
            ds = ds + as
            !omp end parallel workshare

            call field_decompose_semi_spectral(ds)
!            call fftxys2p(ds, df)
            
          end subroutine spectral_diffz

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a divergent flow field (ud, vd, wd) = grad(phi) where
        ! Lap(phi) = div (given).
        subroutine diverge(div,  ud, vd, wd)
            double precision, intent(inout)  :: div(0:nz, ny, nx)
            double precision, intent(out)    :: ud(0:nz, ny, nx), vd(0:nz, ny, nx), wd(0:nz, ny, nx)
            double precision                 :: ds(0:nz, nx, ny)
            double precision                 :: us(0:nz, nx, ny), vs(0:nz, nx, ny), ws(0:nz, nx, ny)

            !------------------------------------------------------------------
            ! Convert phi to spectral space (in x & y) as ds:
            call fftxyp2s(div, ds)

            ds(:, 1, 1) = zero

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
            call diffz(ds, ws)

            ! Set vertical boundary values to zero
            ws(0,  :, :) = zero
            ws(nz, :, :) = zero

            ! Reverse FFT to define z velocity component wd:
            call fftxys2p(ws, wd)

        end subroutine

end module inversion_mod
