module rk4_utils
    use parcel_ellipsoid, only : get_B33
    use fields, only : velgradg, tbuoyg, vortg
    use constants, only : zero, one, two, f12
    use parameters, only : nx, ny, nz, dxi
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
            double precision, intent(in) :: Bin(5)
            double precision, intent(in) :: S(5)
            double precision, intent(in) :: vorticity(3)
            double precision, intent(in) :: volume
            double precision             :: Bout(5), B33
            double precision             :: dudx, dudy, dudz, &
                                            dvdx, dvdy, dvdz, &
                                            dwdx, dwdy, dwdz

            dudx = S(1)
            dudy = S(2)
            dvdy = S(3)
            dwdx = S(4)
            dwdy = S(5)

            ! du/dz = \omegay + dw/dx
            dudz = vorticity(2) + dwdx

            ! dv/dx \omegaz + du/dy
            dvdx = vorticity(3) + dudy

            ! dv/dz = dw/dy - \omegax
            dvdz = dvdy - vorticity(1)

            ! dw/dz = - (du/dx + dv/dy)
            dwdz = - (dudx + dvdy)

            B33 = get_B33(Bin, volume)

            ! dB11/dt = 2 * (du/dx * B11 + du/dy * B12 + du/dz * B13)
            Bout(1) = two * (dudx * Bin(1) + dudy * Bin(2) + dudz * Bin(3))

            ! dB12/dt =
            Bout(2) = dvdx * Bin(1) & !   dv/dx * B11
                    - dwdz * Bin(2) & ! - dw/dz * B12
                    + dvdz * Bin(3) & ! + dv/dz * B13
                    + dudy * Bin(4) & ! + du/dy * B22
                    + dudz * Bin(5)   ! + du/dz * B23

            ! dB13/dt =
            Bout(3) = dwdx * Bin(1) & !   dw/dx * B11
                    + dwdy * Bin(2) & ! + dw/dy * B12
                    - dvdy * Bin(3) & ! - dv/dy * B13
                    + dudy * Bin(5) & ! + du/dy * B23
                    + dudz * B33      ! + du/dz * B33

            ! dB22/dt = 2 * (dv/dx * B12 + dv/dy * B22 + dv/dz * B23)
            Bout(4) = two * (dvdx * Bin(2) + dvdy * Bin(4) + dvdz * Bin(5))

            ! dB23/dt =
            Bout(5) = dwdx * Bin(2) & !   dw/dx * B12
                    + dvdx * Bin(3) & ! + dv/dx * B13
                    + dwdy * Bin(4) & ! + dw/dy * B22
                    - dudx * Bin(5) & ! - du/dx * B23
                    + dvdz * B33      ! + dv/dz * B33
        end function get_dBdt

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax, strain(3, 3), D(3)
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
                        !    dv/dx = \omegaz + du/dy
                        !    du/dz = \omegay + dw/dx
                        !    dv/dz = dw/dy - \omegax
                        !    dw/dz = - (du/dx + dv/dy)
                        strain(1, 1) = velgradg(iz, iy, ix, 1)                              ! du/dx
                        strain(1, 2) = two * velgradg(iz, iy, ix, 2) + vortg(iz, iy, ix, 3) ! du/dy + dv/dx
                        strain(1, 3) = two * velgradg(iz, iy, ix, 4) + vortg(iz, iy, ix, 2) ! du/dz + dw/dx
                        strain(2, 2) = velgradg(iz, iy, ix, 3)                              ! dv/dy
                        strain(2, 3) = two * velgradg(iz, iy, ix, 5) - vortg(iz, iy, ix, 1) ! dv/dz + dw/dy
                        strain(3, 3) = -(velgradg(iz, iy, ix, 1) + velgradg(iz, iy, ix, 3)) ! dw/dz

                        ! calculate its eigenvalues (strain is overwritten and diagonal entries
                        ! will be the eigenvalues sorted in descending order), i.e.
                        ! the largest eigenvalue is in D(1). The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call jacobi_eigenvalues(strain, D)

                        gmax = max(gmax, D(1))
                    enddo
                enddo
            enddo

            !
            ! buoyancy gradient (central difference)
            !

            ! db/dx inner part
            gradb(:, :, 1:nx-2) = f12 * dxi(1) * (tbuoyg(0:nz, :, 2:nx-1) - tbuoyg(0:nz, :, 0:nx-3))

            ! db/dx boundary grid points (make use of periodicity)
            gradb(:, :, 0)    = f12 * dxi(1) * (tbuoyg(0:nz, :, 1) - tbuoyg(0:nz, :, nx-1))
            gradb(:, :, nx-1) = f12 * dxi(1) * (tbuoyg(0:nz, :, 0) - tbuoyg(0:nz, :, nx-2))

            db2 = gradb ** 2

            ! db/dy inner part
            gradb(:, 1:ny-2, :) = f12 * dxi(2) * (tbuoyg(0:nz, 2:ny-1, :) - tbuoyg(0:nz, 0:ny-3, :))

            ! db/dy boundary grid points (make use of periodicity)
            gradb(:, 0, :)    = f12 * dxi(2) * (tbuoyg(0:nz, 1, :) - tbuoyg(0:nz, ny-1, :))
            gradb(:, ny-1, :) = f12 * dxi(2) * (tbuoyg(0:nz, 0, :) - tbuoyg(0:nz, ny-2, :))

            db2 = db2 + gradb ** 2

            ! db/dz
            gradb = f12 * dxi(3) * (tbuoyg(1:nz+1, :, :) - tbuoyg(-1:nz-1, :, :))

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
