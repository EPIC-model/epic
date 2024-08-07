module rk4_utils
    use parcel_ellipse, only : get_B22
    use fields, only : velgradg, tbuoyg, vtend
    use constants, only : zero, one, two, f12, f14, f23, c_matexp_x1, c_matexp_x2, c_matexp_x4, &
                          c_matexp_x5, c_matexp_x6, c_matexp_x7, c_matexp_y2
    use parameters, only : nx, nz, dxi
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif

    implicit none
    double precision, parameter :: Imat(2,2)=reshape((/one,zero,zero,one/), (/2,2/))

    contains

        subroutine evolve_ellipsoid(B, S, volume, dt_sub)
            double precision, intent(inout) :: B(2)
            double precision, intent(in) :: S(4)
            double precision, intent(in) :: volume
            double precision, intent(in) :: dt_sub
            double precision :: Bmat(2,2)
            double precision :: Smat(2,2)
            double precision :: Qmat(2,2)
            double precision :: Rmat(2,2)
            double precision :: Rmat2(2,2)
            double precision :: Rmat4(2,2)
            double precision :: Rmat8(2,2)

            Bmat(1, 1) = B(1) ! B11
            Bmat(1, 2) = B(2) ! B12
            Bmat(2, 1) = B(1) ! B21
            Bmat(2, 2) = get_B22(B(1), B(2), volume)

            Smat(1, 1) = S(1) ! S11
            Smat(1, 2) = S(2) ! S12
            Smat(2, 1) = S(3) ! S21
            Smat(2, 2) = S(4) ! S33

            ! Bader, P., Blanes, S., & Casas, F. (2019). 
            ! Computing the matrix exponential with an optimized Taylor polynomial approximation. 
            ! Mathematics, 7(12), 1174.
            ! Using 8th order Taylor with 2 ward steps
            ! Possibly this is overkill and ward steps can be removed
            ! This does not save much computation though

            ! Bader, P., Blanes, S., & Casas, F. (2019). 
            ! Computing the matrix exponential with an optimized Taylor polynomial approximation. 
            ! Mathematics, 7(12), 1174.
            ! Using 8th order Taylor with 2 ward steps
            ! Possibly this is overkill and ward steps can be removed
            ! This does not save much computation though

            Rmat = f14 * dt_sub * Smat
            Rmat2 = matmul(Rmat, Rmat)
            Rmat4 = matmul(Rmat2, c_matexp_x1 * Rmat + c_matexp_x2 * Rmat2)
            Rmat8 = matmul(f23 * Rmat2 + Rmat4, c_matexp_x4 * Imat + c_matexp_x5 * Rmat + c_matexp_x6 * Rmat2 + c_matexp_x7 * Rmat4)
            Qmat = Imat + Rmat + c_matexp_y2 * Rmat2 + Rmat8
            Qmat = matmul(Qmat, Qmat)
            Qmat = matmul(Qmat, Qmat)
            Bmat = matmul(Qmat, matmul(Bmat, transpose(Qmat)))

            B(1) = Bmat(1, 1)
            B(2) = Bmat(1, 2)

        end subroutine evolve_ellipsoid

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
            gmax = f12 * dsqrt(maxval((velgradg(0:nz, :, 1) - velgradg(0:nz, :, 4)) ** 2 + &
                                        (velgradg(0:nz, :, 2) + velgradg(0:nz, :, 3)) ** 2))
            gmax = max(epsilon(gmax), gmax)

            ! buoyancy gradient

            ! db/dz (central difference)
            dbdz(0:nz, :) = f12 * dxi(2) * (tbuoyg(1:nz+1, :) - tbuoyg(-1:nz-1, :))

            bmax = dsqrt(dsqrt(maxval(vtend(0:nz, :) ** 2 + dbdz ** 2)))
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
