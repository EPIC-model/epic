module rk4_utils
    use parcel_ellipse, only : get_B22
    use fields, only : velgradg, tbuoyg, vtend
    use constants, only : zero, one, two, f12
    use parameters, only : nx, nz, dxi

    implicit none

    contains

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


        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax
            double precision             :: dbdz(0:nz, 0:nx-1)

            if (time%is_adaptive) then
                ! velocity strain
                gmax = f12 * dsqrt(maxval((velgradg(0:nz, :, 1) - velgradg(0:nz, :, 4)) ** 2 + &
                                          (velgradg(0:nz, :, 2) + velgradg(0:nz, :, 3)) ** 2))
                gmax = max(epsilon(gmax), gmax)

                ! buoyancy gradient

                ! db/dz (central difference)
                dbdz(0:nz, :) = f12 * dxi(2) * (tbuoyg(1:nz+1, :) - tbuoyg(-1:nz-1, :))

                bmax = dsqrt(dsqrt(maxval(vtend(0:nz, :) ** 2 + dbdz ** 2)))
                bmax = max(epsilon(bmax), bmax)

                dt = min(time%alpha_s / gmax, time%alpha_b / bmax)
            else
                dt = time%dt
            endif

            if (time%precise_stop .and. (t + dt > time%limit)) then
                dt = time%limit - t
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step

end module rk4_utils
