module rk4_utils
    use parcel_ellipsoid, only : get_B33
    use fields, only : velgradg, tbuoyg, vortg
    use constants, only : zero, one, two, f12
    use parameters, only : nx, ny, nz, dxi
    use jacobi, only : jacobi_diagonalise
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif

    implicit none

    contains

        ! Advance the B matrix.
        ! @param[in] Bin are the B matrix components of the parcel
        ! @param[in] S is the local velocity strain
        ! @param[in] volume is the parcel volume
        ! @returns dB/dt in Bout
        function get_dBdt(Bin, S, volume) result(Bout)
            double precision, intent(in) :: Bin(5)
            double precision, intent(in) :: S(9)
            double precision, intent(in) :: volume
            double precision             :: Bout(5), B33

            ! du/dx = S(1)
            ! du/dy = S(2)
            ! du/dz = S(3)
            ! dv/dx = S(4)
            ! dv/dy = S(5)
            ! dv/dz = S(6)
            ! dw/dx = S(7)
            ! dw/dy = S(8)
            ! dw/dz = S(9)

            B33 = get_B33(Bin, volume)

            ! dB11/dt = 2 * (du/dx * B11 + du/dy * B12 + du/dz * B13)
            Bout(1) = two * (S(1) * Bin(1) + S(2) * Bin(2) + S(3) * Bin(3))

            ! dB12/dt =
            Bout(2) = S(4) * Bin(1) & !   dv/dx * B11
                    - S(9) * Bin(2) & ! - dw/dz * B12
                    + S(6) * Bin(3) & ! + dv/dz * B13
                    + S(2) * Bin(4) & ! + du/dy * B22
                    + S(3) * Bin(5)   ! + du/dz * B23

            ! dB13/dt =
            Bout(3) = S(7) * Bin(1) & !   dw/dx * B11
                    + S(8) * Bin(2) & ! + dw/dy * B12
                    - S(5) * Bin(3) & ! - dv/dy * B13
                    + S(2) * Bin(5) & ! + du/dy * B23
                    + S(3) * B33      ! + du/dz * B33

            ! dB22/dt = 2 * (dv/dx * B12 + dv/dy * B22 + dv/dz * B23)
            Bout(4) = two * (S(4) * Bin(2) + S(5) * Bin(4) + S(6) * Bin(5))

            ! dB23/dt =
            Bout(5) = S(7) * Bin(2) & !   dw/dx * B12
                    + S(4) * Bin(3) & ! + dv/dx * B13
                    + S(8) * Bin(4) & ! + dw/dy * B22
                    - S(1) * Bin(5) & ! - du/dx * B23
                    + S(6) * B33      ! + dv/dz * B33
        end function get_dBdt

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax, strain(3, 3)
            double precision             :: dbdx(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dbdy(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dbdz(0:nz, 0:ny-1, 0:nx-1)
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

                        strain(2, 1) = strain(1, 2)
                        strain(3, 1) = strain(1, 3)
                        strain(3, 2) = strain(2, 3)

                        ! calculate its eigenvalues (strain is overwritten and will be
                        ! the diagonal matrix with the eigenvalues sorted in descending order), i.e.
                        ! the largest eigenvalue is in strain(1, 1).
                        call jacobi_diagonalise(strain)

                        gmax = max(gmax, strain(1, 1))
                    enddo
                enddo
            enddo

            !
            ! buoyancy gradient (central difference)
            !

            ! db/dx inner part
            dbdx(:, :, 1:nx-2) = f12 * dxi(1) * (tbuoyg(0:nz, :, 2:nx-1) - tbuoyg(0:nz, :, 0:nx-3))

            ! db/dx boundary grid points (make use of periodicity)
            dbdx(:, :, 0)    = f12 * dxi(1) * (tbuoyg(0:nz, :, 1) - tbuoyg(0:nz, :, nx-1))
            dbdx(:, :, nx-1) = f12 * dxi(1) * (tbuoyg(0:nz, :, 0) - tbuoyg(0:nz, :, nx-2))

            ! db/dy inner part
            dbdy(:, 1:ny-2, :) = f12 * dxi(2) * (tbuoyg(0:nz, 2:ny-1, :) - tbuoyg(0:nz, 0:ny-3, :))

            ! db/dy boundary grid points (make use of periodicity)
            dbdy(:, 0, :)    = f12 * dxi(2) * (tbuoyg(0:nz, 1, :) - tbuoyg(0:nz, ny-1, :))
            dbdy(:, ny-1, :) = f12 * dxi(2) * (tbuoyg(0:nz, 0, :) - tbuoyg(0:nz, ny-2, :))

            ! db/dz
            dbdz = f12 * dxi(3) * (tbuoyg(1:nz+1, :, :) - tbuoyg(-1:nz-1, :, :))


            bmax = dsqrt(dsqrt(maxval(dbdx ** 2 + dbdy ** 2 + dbdz ** 2)))
            bmax = max(epsilon(bmax), bmax)

            dt = min(time%alpha / gmax, time%alpha / bmax)
#ifdef ENABLE_VERBOSE
            fname = trim(output%h5_basename) // '_alpha_time_step.asc'
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
