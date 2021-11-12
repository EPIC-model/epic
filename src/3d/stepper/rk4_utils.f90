module rk4_utils
    use parcel_ellipsoid, only : get_B33
    use fields, only : velgradg, tbuoyg, vtend
    use constants, only : zero, one, two, f12
    use parameters, only : nx, nz, dxi
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif

    implicit none

    contains

        ! Advance the B matrix.
        ! @param[in] Bin are the B matrix components of the parcel
        ! @param[in] S is the local velocity strain
        ! @param[in] volume is the parcel volume
        ! @returns the updated B matrix components (B11 and B12) in Bout
        function get_B(Bin, S, volume) result(Bout)
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
            Bout(4) = two * (S(4) * Bin(3) + S(5) * Bin(4) + S(6) * Bin(5))

            ! dB23/dt =
            Bout(5) = S(7) * Bin(2) & !   dw/dx * B12
                    + S(4) * Bin(3) & ! + dv/dx * B13
                    + S(8) * Bin(4) & ! + dw/dy * B22
                    - S(1) * Bin(5) & ! - du/dx * B23
                    + S(6) * B33      ! + dv/dz * B33
        end function get_B

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax
            double precision             :: dbdz(0:nz, 0:nx-1)
#if ENABLE_VERBOSE
            logical                      :: exists = .false.
            character(:), allocatable    :: fname
#endif

            ! velocity strain
            gmax = 1.0 !FIXME f12 * dsqrt(maxval((velgradg(0:nz, :, 1) - velgradg(0:nz, :, 4)) ** 2 + &
!                                         (velgradg(0:nz, :, 2) + velgradg(0:nz, :, 3)) ** 2))
            gmax = max(epsilon(gmax), gmax)

            ! buoyancy gradient

            ! db/dz (central difference)
            dbdz(0:nz, :) = 1.0 !FIXME f12 * dxi(2) * (tbuoyg(1:nz+1, :) - tbuoyg(-1:nz-1, :))

            bmax = 1.0 !FIXME dsqrt(dsqrt(maxval(vtend(0:nz, :) ** 2 + dbdz ** 2)))
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
