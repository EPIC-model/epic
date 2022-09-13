module rk4_utils
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use parcel_container, only : I_B11, I_B12, I_B13, I_B22, I_B23
    use parcel_ellipsoid, only : get_B33
    use fields, only : velgradg, tbuoyg, vortg, I_DUDX, I_DUDY, I_DVDY, I_DWDX, I_DWDY
    use constants, only : zero, one, two, f12
    use parameters, only : nx, ny, nz, dxi, vcell
    use jacobi, only : jacobi_eigenvalues
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif

    implicit none

    contains

        ! Advance the B matrix.
        ! @param[in] Bin are the B matrix components of the parcel
        ! @param[in] S is the local velocity strain
        ! @param[in] vorticity of parcel
        ! @param[in] volume is the parcel volume
        ! @returns dB/dt in Bout
        function get_dBdt(Bin, S, vorticity, volume) result(Bout)
            double precision, intent(in) :: Bin(I_B23)
            double precision, intent(in) :: S(5)
            double precision, intent(in) :: vorticity(n_dim)
            double precision, intent(in) :: volume
            double precision             :: Bout(5), B33
            double precision             :: dudz, dvdx, dvdz, dwdz

            ! du/dz = \eta + dw/dx
            dudz = vorticity(I_Y) + S(I_DWDX)

            ! dv/dx \zeta + du/dy
            dvdx = vorticity(I_Z) + S(I_DUDY)

            ! dv/dz = dw/dy - \xi
            dvdz = S(I_DWDY) - vorticity(I_X)

            ! dw/dz = - (du/dx + dv/dy)
            dwdz = - (S(I_DUDX) + S(I_DVDY))

            B33 = get_B33(Bin, volume)

            ! dB11/dt = 2 * (du/dx * B11 + du/dy * B12 + du/dz * B13)
            Bout(I_B11) = two * (S(I_DUDX) * Bin(I_B11) + S(I_DUDY) * Bin(I_B12) + dudz * Bin(I_B13))

            ! dB12/dt =
            Bout(I_B12) = dvdx      * Bin(I_B11) & !   dv/dx * B11
                        - dwdz      * Bin(I_B12) & ! - dw/dz * B12
                        + dvdz      * Bin(I_B13) & ! + dv/dz * B13
                        + S(I_DUDY) * Bin(I_B22) & ! + du/dy * B22
                        + dudz      * Bin(I_B23)   ! + du/dz * B23

            ! dB13/dt =
            Bout(I_B13) = S(I_DWDX) * Bin(I_B11) & !   dw/dx * B11
                        + S(I_DWDY) * Bin(I_B12) & ! + dw/dy * B12
                        - S(I_DVDY) * Bin(I_B13) & ! - dv/dy * B13
                        + S(I_DUDY) * Bin(I_B23) & ! + du/dy * B23
                        + dudz      * B33         ! + du/dz * B33

            ! dB22/dt = 2 * (dv/dx * B12 + dv/dy * B22 + dv/dz * B23)
            Bout(I_B22) = two * (dvdx * Bin(I_B12) + S(I_DVDY) * Bin(I_B22) + dvdz * Bin(I_B23))

            ! dB23/dt =
            Bout(I_B23) = S(I_DWDX) * Bin(I_B12) & !   dw/dx * B12
                        + dvdx      * Bin(I_B13) & ! + dv/dx * B13
                        + S(I_DWDY) * Bin(I_B22) & ! + dw/dy * B22
                        - S(I_DUDX) * Bin(I_B23) & ! - du/dx * B23
                        + dvdz      * B33          ! + dv/dz * B33
        end function get_dBdt

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax, strain(n_dim, n_dim), D(n_dim)
            double precision             :: gradb(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: db2(0:nz, 0:ny-1, 0:nx-1)
            integer                      :: ix, iy, iz
#if ENABLE_VERBOSE
            logical                      :: exists = .false.
            character(:), allocatable    :: fname
#endif

            !
            ! velocity strain
            !
            ! find largest stretch -- this corresponds to largest
            ! eigenvalue over all local symmetrised strain matrices.
            gmax = epsilon(gmax)
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 0, nz
                        ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
                        ! where
                        !     /u_x u_y u_z\
                        ! S = |v_x v_y v_z|
                        !     \w_x w_y w_z/
                        ! with u_* = du/d* (also derivatives of v and w).
                        ! The derivatives dv/dx, du/dz, dv/dz and dw/dz are calculated
                        ! with vorticity or the assumption of incompressibility
                        ! (du/dx + dv/dy + dw/dz = 0):
                        !    dv/dx = \zeta + du/dy
                        !    du/dz = \eta + dw/dx
                        !    dv/dz = dw/dy - \xi
                        !    dw/dz = - (du/dx + dv/dy)
                        !
                        !                         /  2 * u_x  u_y + v_x u_z + w_x\
                        ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
                        !                         \u_z + w_x  v_z + w_y   2 * w_z/
                        !
                        ! S11 = du/dx
                        ! S12 = 1/2 * (du/dy + dv/dx) = 1/2 * (2 * du/dy + \zeta) = du/dy + 1/2 * \zeta
                        ! S13 = 1/2 * (du/dz + dw/dx) = 1/2 * (\eta + 2 * dw/dx) = 1/2 * \eta + dw/dx
                        ! S22 = dv/dy
                        ! S23 = 1/2 * (dv/dz + dw/dy) = 1/2 * (2 * dw/dy - \xi) = dw/dy - 1/2 * \xi
                        ! S33 = dw/dz = - (du/dx + dv/dy)
                        strain(1, 1) = velgradg(iz, iy, ix, I_DUDX)                                   ! S11
                        strain(1, 2) = velgradg(iz, iy, ix, I_DUDY) + f12 * vortg(iz, iy, ix, I_Z)    ! S12
                        strain(1, 3) = velgradg(iz, iy, ix, I_DWDX) + f12 * vortg(iz, iy, ix, I_Y)    ! S13
                        strain(2, 2) = velgradg(iz, iy, ix, I_DVDY)                                   ! S22
                        strain(2, 3) = velgradg(iz, iy, ix, I_DWDY) - f12 * vortg(iz, iy, ix, I_X)    ! S23
                        strain(3, 3) = -(velgradg(iz, iy, ix, I_DUDX) + velgradg(iz, iy, ix, I_DVDY)) ! S33

                        ! calculate its eigenvalues. The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call jacobi_eigenvalues(strain, D)

                        ! we must take the largest eigenvalue in magnitude (absolute value)
                        gmax = max(gmax, maxval(abs(D)))
                    enddo
                enddo
            enddo

            !
            ! buoyancy gradient (central difference)
            !

            ! db/dx inner part
            gradb(:, :, 1:nx-2) = f12 * dxi(I_X) * (tbuoyg(0:nz, :, 2:nx-1) - tbuoyg(0:nz, :, 0:nx-3))

            ! db/dx boundary grid points (make use of periodicity)
            gradb(:, :, 0)    = f12 * dxi(I_X) * (tbuoyg(0:nz, :, 1) - tbuoyg(0:nz, :, nx-1))
            gradb(:, :, nx-1) = f12 * dxi(I_X) * (tbuoyg(0:nz, :, 0) - tbuoyg(0:nz, :, nx-2))

            db2 = gradb ** 2

            ! db/dy inner part
            gradb(:, 1:ny-2, :) = f12 * dxi(I_Y) * (tbuoyg(0:nz, 2:ny-1, :) - tbuoyg(0:nz, 0:ny-3, :))

            ! db/dy boundary grid points (make use of periodicity)
            gradb(:, 0, :)    = f12 * dxi(I_Y) * (tbuoyg(0:nz, 1, :) - tbuoyg(0:nz, ny-1, :))
            gradb(:, ny-1, :) = f12 * dxi(I_Y) * (tbuoyg(0:nz, 0, :) - tbuoyg(0:nz, ny-2, :))

            db2 = db2 + gradb ** 2

            ! db/dz
            gradb = f12 * dxi(I_Z) * (tbuoyg(1:nz+1, :, :) - tbuoyg(-1:nz-1, :, :))

            bmax = dsqrt(dsqrt(maxval(db2 + gradb ** 2)))
            bmax = max(epsilon(bmax), bmax)

            dt = min(time%alpha / gmax, time%alpha / bmax)
#ifdef ENABLE_VERBOSE
            fname = trim(output%basename) // '_alpha_time_step.asc'
            inquire(file=fname, exist=exists)
            ! 23 August
            ! https://stackoverflow.com/questions/15526203/single-command-to-open-a-file-or-create-it-and-the-append-data
            if ((t /= zero) .and. exists) then
                open(unit=1235, file=fname, status='old', position='append')
            else
                open(unit=1235, file=fname, status='replace')
                write(1235, *) '  # time (s)                \alpha_s/\gamma_{max}     \alpha_b/N_{max}'
            endif

            write(1235, *) t, time%alpha / gmax, time%alpha / bmax

            close(1235)
#endif
            if (time%precise_stop .and. (t + dt > time%limit)) then
                dt = time%limit - t
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step

end module rk4_utils
