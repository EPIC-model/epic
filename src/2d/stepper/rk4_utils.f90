module rk4_utils
    use parcel_ellipse, only : get_B22
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
            double precision, intent(in) :: Bin(2)
            double precision, intent(in) :: S(4)
            double precision, intent(in) :: volume
            double precision             :: Bout(2)

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(1) = two * (S(1) * Bin(1) + S(2) * Bin(2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(2) = S(3) * Bin(1) + S(2) * get_B22(Bin(1), Bin(2), volume)

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
            character(512)               :: fname
#endif

            ! velocity strain
            gmax = f12 * sqrt(maxval((velgradg(0:nz, :, 1) - velgradg(0:nz, :, 4)) ** 2 + &
                                        (velgradg(0:nz, :, 2) + velgradg(0:nz, :, 3)) ** 2))
            gmax = max(epsilon(gmax), gmax)

            ! buoyancy gradient

            ! db/dz (central difference)
            dbdz(0:nz, :) = f12 * dxi(2) * (tbuoyg(1:nz+1, :) - tbuoyg(-1:nz-1, :))

            bmax = sqrt(sqrt(maxval(vtend(0:nz, :) ** 2 + dbdz ** 2)))
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
