module rk_utils
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use parcel_ellipsoid, only : I_B11, I_B12, I_B13, I_B22, I_B23
    use fields, only : velgradg, tbuoyg, vortg, I_DUDX, I_DUDY, I_DVDY, I_DWDX, I_DWDY, strain_mag
    use field_mpi, only : field_halo_fill_scalar
    use constants, only : zero, one, two, f12
    use parameters, only : nx, ny, nz, dxi, vcell
    use scherzinger, only : scherzinger_eigenvalues
    use mpi_layout, only : box
    use mpi_environment
    use mpi_utils, only : mpi_exit_on_error
    use parcels_mod, only : parcels
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif
#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
    use physics, only : bfsq
#endif

    implicit none

    contains

        ! Advance the B matrix.
        ! @param[in] n parcel index
        ! @returns dB/dt in Bout
        function get_dBdt(n) result(Bout)
            integer,         intent(in) :: n
            double precision            :: Bout(5), B33
            double precision            :: dudz, dvdx, dvdz, dwdz

            ! du/dz = \eta + dw/dx
            dudz = parcels%vorticity(I_Y, n) + parcels%strain(I_DWDX, n)

            ! dv/dx \zeta + du/dy
            dvdx = parcels%vorticity(I_Z, n) + parcels%strain(I_DUDY, n)

            ! dv/dz = dw/dy - \xi
            dvdz = parcels%strain(I_DWDY, n) - parcels%vorticity(I_X, n)

            ! dw/dz = - (du/dx + dv/dy)
            dwdz = - (parcels%strain(I_DUDX, n) + parcels%strain(I_DVDY, n))

            B33 = parcels%get_B33(n)

            ! dB11/dt = 2 * (du/dx * B11 + du/dy * B12 + du/dz * B13)
            Bout(I_B11) = two * (parcels%strain(I_DUDX, n) * parcels%B(I_B11, n)  + &
                                 parcels%strain(I_DUDY, n) * parcels%B(I_B12, n)  + &
                                 dudz                      * parcels%B(I_B13, n))

            ! dB12/dt =
            Bout(I_B12) = dvdx                      * parcels%B(I_B11, n) & !   dv/dx * B11
                        - dwdz                      * parcels%B(I_B12, n) & ! - dw/dz * B12
                        + dvdz                      * parcels%B(I_B13, n) & ! + dv/dz * B13
                        + parcels%strain(I_DUDY, n) * parcels%B(I_B22, n) & ! + du/dy * B22
                        + dudz                      * parcels%B(I_B23, n)   ! + du/dz * B23

            ! dB13/dt =
            Bout(I_B13) = parcels%strain(I_DWDX, n) * parcels%B(I_B11, n) & !   dw/dx * B11
                        + parcels%strain(I_DWDY, n) * parcels%B(I_B12, n) & ! + dw/dy * B12
                        - parcels%strain(I_DVDY, n) * parcels%B(I_B13, n) & ! - dv/dy * B13
                        + parcels%strain(I_DUDY, n) * parcels%B(I_B23, n) & ! + du/dy * B23
                        + dudz                      * B33                   ! + du/dz * B33

            ! dB22/dt = 2 * (dv/dx * B12 + dv/dy * B22 + dv/dz * B23)
            Bout(I_B22) = two * (dvdx                      * parcels%B(I_B12, n) + &
                                 parcels%strain(I_DVDY, n) * parcels%B(I_B22, n) + &
                                 dvdz                      * parcels%B(I_B23, n))

            ! dB23/dt =
            Bout(I_B23) = parcels%strain(I_DWDX, n) * parcels%B(I_B12, n) & !   dw/dx * B12
                        + dvdx                      * parcels%B(I_B13, n) & ! + dv/dx * B13
                        + parcels%strain(I_DWDY, n) * parcels%B(I_B22, n) & ! + dw/dy * B22
                        - parcels%strain(I_DUDX, n) * parcels%B(I_B23, n) & ! - du/dx * B23
                        + dvdz                      * B33                   ! + dv/dz * B33
        end function get_dBdt

        ! Calculate velocity strain
        ! @param[in] velocity gradient tensor at grid point
        ! @param[in] vorticity at grid point
        ! @returns 3x3 strain matrix
        function get_strain(velgradgp, vortgp) result(strain)
            double precision, intent(in) :: velgradgp(5)
            double precision, intent(in) :: vortgp(n_dim)
            double precision             :: strain(3,3)

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

            strain(1, 1) = velgradgp(I_DUDX)                        ! S11
            strain(1, 2) = velgradgp(I_DUDY) + f12 * vortgp(I_Z)    ! S12
            strain(1, 3) = velgradgp(I_DWDX) + f12 * vortgp(I_Y)    ! S13
            strain(2, 1) = strain(1, 2)
            strain(2, 2) = velgradgp(I_DVDY)                        ! S22
            strain(2, 3) = velgradgp(I_DWDY) - f12 * vortgp(I_X)    ! S23
            strain(3, 1) = strain(1, 3)
            strain(3, 2) = strain(2, 3)
            strain(3, 3) = -(velgradgp(I_DUDX) + velgradgp(I_DVDY)) ! S33
        end function get_strain

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax, strain(n_dim, n_dim), D(n_dim), local_max(2)
            double precision             :: gradb(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision             :: db2(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz
#if ENABLE_VERBOSE
            logical                      :: exists = .false.
            character(512)               :: fname
#endif

            !
            ! velocity strain
            !
            ! find largest stretch -- this corresponds to largest
            ! eigenvalue over all local symmetrised strain matrices.
            gmax = epsilon(gmax)
            do ix = box%lo(1), box%hi(1)
                do iy = box%lo(2), box%hi(2)
                    do iz = 0, nz
                        strain = get_strain(velgradg(iz, iy, ix,:), vortg(iz, iy, ix, :))
                        ! calculate its eigenvalues. The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call scherzinger_eigenvalues(strain, D)

                        ! we must take the largest eigenvalue in magnitude (absolute value)
                        gmax = max(gmax, maxval(abs(D)))
                    enddo
                enddo
            enddo

            !
            ! buoyancy gradient (central difference)
            !

            ! ensure that halo grid points are filled
            call field_halo_fill_scalar(tbuoyg, l_alloc=.true.)

            ! db/dx
            gradb = f12 * dxi(I_X) * (tbuoyg(0:nz, box%lo(2):box%hi(2), box%lo(1)+1:box%hi(1)+1) &
                                    - tbuoyg(0:nz, box%lo(2):box%hi(2), box%lo(1)-1:box%hi(1)-1))

            db2 = gradb ** 2

            ! db/dy
            gradb = f12 * dxi(I_Y) * (tbuoyg(0:nz, box%lo(2)+1:box%hi(2)+1, box%lo(1):box%hi(1)) &
                                    - tbuoyg(0:nz, box%lo(2)-1:box%hi(2)-1, box%lo(1):box%hi(1)))

            db2 = db2 + gradb ** 2

            ! db/dz
            gradb = f12 * dxi(I_Z) * (tbuoyg(1:nz+1,  box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                                    - tbuoyg(-1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
            gradb = gradb + bfsq
#endif

            bmax = dsqrt(dsqrt(maxval(db2 + gradb ** 2)))
            bmax = max(epsilon(bmax), bmax)

            local_max(1) = gmax
            local_max(2) = bmax

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               local_max(1:2),          &
                               2,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            gmax = local_max(1)
            bmax = local_max(2)

            dt = min(time%alpha / gmax, time%alpha / bmax)
#ifdef ENABLE_VERBOSE
            if (world%rank == world%root) then
                fname = trim(output%basename) // '_alpha_time_step.asc'
                inquire(file=trim(fname), exist=exists)
                ! 23 August
                ! https://stackoverflow.com/questions/15526203/single-command-to-open-a-file-or-create-it-and-the-append-data
                if ((t /= zero) .and. exists) then
                    open(unit=1235, file=trim(fname), status='old', position='append')
                else
                    open(unit=1235, file=trim(fname), status='replace')
                    write(1235, *) '  # time (s)                \alpha_s/\gamma_{max}     \alpha_b/N_{max}'
                endif

                write(1235, *) t, time%alpha / gmax, time%alpha / bmax

                close(1235)
            endif
#endif
            if (time%precise_stop .and. (t + dt > time%limit)) then
                dt = time%limit - t
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                call mpi_exit_on_error
            endif
        end function get_time_step

        ! @returns the strain magnitude for use in damping routine
        subroutine get_strain_magnitude_field
            double precision             :: strain(n_dim, n_dim)
            integer                      :: ix, iy, iz

            do ix = box%lo(1), box%hi(1)
                do iy = box%lo(2), box%hi(2)
                    do iz = 0, nz
                        strain = get_strain(velgradg(iz, iy, ix,:), vortg(iz, iy, ix, :))
                        strain_mag(iz, iy, ix) = sqrt(two * (strain(1, 1) * strain(1, 1) +&
                                                             strain(1, 2) * strain(1, 2) +&
                                                             strain(1, 3) * strain(1, 3) +&
                                                             strain(2, 1) * strain(2, 1) +&
                                                             strain(2, 2) * strain(2, 2) +&
                                                             strain(2, 3) * strain(2, 3) +&
                                                             strain(3, 1) * strain(3, 1) +&
                                                             strain(3, 2) * strain(3, 2) +&
                                                             strain(3, 3) * strain(3, 3) ))
                   enddo
                   ! Reflect beyond boundaries to ensure damping is conservative
                   ! This is because the points below the surface contribute to the level above
                   strain_mag(-1, iy, ix) = strain_mag(1, iy, ix)
                   strain_mag(nz+1, iy, ix) = strain_mag(nz-1, iy, ix)
                enddo
            enddo

          ! We need this halo fill to obtain good conservation
          call field_halo_fill_scalar(strain_mag, l_alloc=.true.)

        end subroutine get_strain_magnitude_field

end module rk_utils
