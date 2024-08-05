module rk_utils
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use parcel_ellipsoid, only : get_B33, I_B11, I_B12, I_B13, I_B22, I_B23, I_B33
    use fields, only : velgradg, tbuoyg, vortg, I_DUDX, I_DUDY, I_DUDZ, I_DVDX, I_DVDY, I_DVDZ, I_DWDX, I_DWDY, strain_mag
    use field_mpi, only : field_halo_fill_scalar
    use constants, only : zero, one, two, f12, f14, f23, xx1, xx2, xx4, xx5, xx6, xx7, yy2
    use parameters, only : nx, ny, nz, dxi, vcell
    use scherzinger, only : scherzinger_eigenvalues
    use mpi_layout, only : box
    use mpi_environment
    use mpi_utils, only : mpi_exit_on_error
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif
#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
    use physics, only : bfsq
#endif

    implicit none
    double precision, parameter :: Imat(3,3) = reshape((/one,zero,zero,zero,one,zero,zero,zero,one/), (/3,3/))

    contains

        subroutine evolve_ellipsoid(B, S, dt_sub)
            double precision, intent(inout) :: B(I_B33)
            double precision, intent(in) :: S(8)
            double precision, intent(in) :: dt_sub
            double precision :: Bmat(3,3)
            double precision :: Smat(3,3)
            double precision :: Rmat(3,3)
            double precision :: Rmat2(3,3)
            double precision :: Rmat4(3,3)
            double precision :: Rmat8(3,3)
            double precision :: Qmat(3,3)

            Bmat(1, 1) = B(I_B11) ! B11
            Bmat(1, 2) = B(I_B12) ! B12
            Bmat(1, 3) = B(I_B13) ! B13
            Bmat(2, 1) = Bmat(1, 2) ! B21
            Bmat(2, 2) = B(I_B22) ! B22
            Bmat(2, 3) = B(I_B23) ! B23
            Bmat(3, 1) = Bmat(1, 3) ! B31
            Bmat(3, 2) = Bmat(2, 3) ! B32
            Bmat(3, 3) = B(I_B33) ! B33

            Smat(1, 1) = S(I_DUDX) ! S11
            Smat(1, 2) = S(I_DUDY) ! S12
            Smat(1, 3) = S(I_DUDZ) ! S13
            Smat(2, 1) = S(I_DVDX) ! S21
            Smat(2, 2) = S(I_DVDY) ! S22
            Smat(2, 3) = S(I_DVDZ) ! S23
            Smat(3, 1) = S(I_DWDX) ! S31
            Smat(3, 2) = S(I_DWDY) ! S32
            Smat(3, 3) = -(Smat(1, 1) + Smat(2, 2)) ! S33

            ! Bader, P., Blanes, S., & Casas, F. (2019). 
            ! Computing the matrix exponential with an optimized Taylor polynomial approximation. 
            ! Mathematics, 7(12), 1174.
            ! Using 8th order Taylor with 2 ward steps
            ! Possibly this is overkill and ward steps can be removed
            ! This does not save much computation though

            Rmat = f14 * dt_sub * Smat
            Rmat2 = matmul(Rmat, Rmat)
            Rmat4 = matmul(Rmat2, xx1 * Rmat + xx2 * Rmat2)
            Rmat8 = matmul(f23 * Rmat2 + Rmat4, xx4 * Imat + xx5 * Rmat + xx6 * Rmat2 + xx7 * Rmat4)
            Qmat = Imat + Rmat + yy2 * Rmat2 + Rmat8
            Qmat = matmul(Qmat, Qmat)
            Qmat = matmul(Qmat, Qmat)
            Bmat = matmul(Qmat, matmul(Bmat, transpose(Qmat)))

            B(I_B11) = Bmat(1, 1)
            B(I_B12) = Bmat(1, 2)
            B(I_B13) = Bmat(1, 3)
            B(I_B22) = Bmat(2, 2)
            B(I_B23) = Bmat(2, 3)
            B(I_B33) = Bmat(3, 3)

        end subroutine


        ! Calculate velocity strain
        ! @param[in] velocity gradient tensor at grid point
        ! @param[in] vorticity at grid point
        ! @returns 3x3 strain matrix
        function get_strain(velgradgp) result(strain)
            double precision, intent(in) :: velgradgp(8)
            double precision             :: strain(3,3)
       
            ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
            ! where
            !     /u_x u_y u_z\
            ! S = |v_x v_y v_z|
            !     \w_x w_y w_z/
            ! with u_* = du/d* (also derivatives of v and w).
            !    dw/dz is obtained from incompressibility
            !    dw/dz = - (du/dx + dv/dy)
            !
            !                         /  2 * u_x  u_y + v_x u_z + w_x\
            ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
            !                         \u_z + w_x  v_z + w_y   2 * w_z/
            !
            ! S11 = du/dx
            ! S12 = 1/2 * (du/dy + dv/dx)
            ! S13 = 1/2 * (du/dz + dw/dx)
            ! S22 = dv/dy
            ! S23 = 1/2 * (dv/dz + dw/dy)
            ! S33 = dw/dz = - (du/dx + dv/dy)

            strain(1, 1) = velgradgp(I_DUDX)                              ! S11
            strain(1, 2) = f12 * (velgradgp(I_DUDY) + velgradgp(I_DVDX))  ! S12
            strain(1, 3) = f12 * (velgradgp(I_DUDZ) + velgradgp(I_DWDX))  ! S13
            strain(2, 1) = strain(1, 2)
            strain(2, 2) = velgradgp(I_DVDY)                              ! S22
            strain(2, 3) = f12 * (velgradgp(I_DVDZ) + velgradgp(I_DWDY))  ! S23
            strain(3, 1) = strain(1, 3)
            strain(3, 2) = strain(2, 3)
            strain(3, 3) = -(velgradgp(I_DUDX) + velgradgp(I_DVDY))       ! S33
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
                        strain = get_strain(velgradg(iz, iy, ix,:))
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
            
            if (time%l_use_fixed_dt) then
                dt = time%fixed_dt
            else
                dt = min(time%alpha / gmax, time%alpha / bmax)
            endif
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
                        strain = get_strain(velgradg(iz, iy, ix,:))
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
