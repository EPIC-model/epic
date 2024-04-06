!===================================================================
!  Uses a divergent flow to push parcels toward or away from grid
!  points so as to bring the area fraction at each grid point back
!  to unity.

!  Here tri-diagonal solves are used in the vertical direction.

!           Initialise data first using inipar.f90
!===================================================================

module surface_parcel_correction
    use sta3dfft, only : fft2d, ifft2d, diff2dx, diff2dy
    use inversion_utils, only : lapinv
    use parcel_ellipse_interpl, only : area2grid
    use interpl, only : bilinear
    use parcels_mod
    use parcel_bc
    use omp_lib
    use mpi_layout, only : box
    use constants
    use parameters, only : acelli, nx, ny, dx, dxi, nz
    use parcel_mpi, only : parcel_communicate
    use field_mpi, only : field_halo_fill_vector
    use mpi_environment
    use mpi_utils, only : mpi_check_for_error
!     use timer, only : start_timer, stop_timer
    use fields, only : volg

    implicit none

    private

!     integer :: surf_lapl_corr_timer, &
!                surf_grad_corr_timer

    ! initial area-weighted vorticity mean
    double precision :: top_vor_bar(3)      &
                      , bot_vor_bar(3)

    public :: apply_surf_laplace,           &
              apply_surf_gradient,          &
              area_correction,             ! &
!               init_surf_parcel_correction,  &
!               apply_surface_vortcor
!               surf_lapl_corr_timer,    &
!               surf_grad_corr_timer


    contains

!         ! Initialise parcel correction module
!         subroutine init_surf_parcel_correction
!             integer          :: n
!             double precision :: bot_sum, top_sum
!             double precision :: buf(8)
!
!             bot_sum = zero
!             top_sum = zero
!             bot_vor_bar = zero
!             top_vor_bar = zero
!
!             !$omp parallel default(shared)
!             !$omp do private(n) reduction(+: bot_vor_bar, bot_sum)
!             do n = 1, bot_parcels%local_num
!                 bot_sum = bot_sum + bot_parcels%area(n)
!                 bot_vor_bar = bot_vor_bar + bot_parcels%vorticity(:, n) * bot_parcels%area(n)
!             enddo
!             !$omp end do
!             !$omp end parallel
!
!             !$omp parallel default(shared)
!             !$omp do private(n) reduction(+: top_vor_bar, top_sum)
!             do n = 1, top_parcels%local_num
!                 top_sum = top_sum + top_parcels%area(n)
!                 top_vor_bar = top_vor_bar + top_parcels%vorticity(:, n) * top_parcels%area(n)
!             enddo
!             !$omp end do
!             !$omp end parallel
!
!             buf(1) = bot_sum
!             buf(2) = top_sum
!             buf(3:5) = bot_vor_bar
!             buf(6:8) = top_vor_bar
!             call MPI_Allreduce(MPI_IN_PLACE,            &
!                                buf(1:8),                &
!                                8,                       &
!                                MPI_DOUBLE_PRECISION,    &
!                                MPI_SUM,                 &
!                                world%comm,              &
!                                world%err)
!
!             call mpi_check_for_error(world, &
!                 "in MPI_Allreduce of parcel_correction::init_surf_parcel_correction.")
!
!             bot_sum = buf(1)
!             top_sum = buf(2)
!             bot_vor_bar = buf(3:5)
!             top_vor_bar = buf(6:8)
!
!             bot_vor_bar = bot_vor_bar / bot_sum
!             top_vor_bar = top_vor_bar / top_sum
!
!         end subroutine init_surf_parcel_correction

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!         subroutine apply_surface_vortcor
!             call m_apply_surface_vortcor(bot_parcels, bot_vor_bar)
!             call m_apply_surface_vortcor(top_parcels, top_vor_bar)
!         end subroutine apply_surface_vortcor

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This subroutine makes sure that the net vorticity
        ! (area-integrated vorticity) remains constant.
        subroutine m_apply_surface_vortcor(pcont, vor_bar)
            type(ellipse_pc_type), intent(inout) :: pcont
            double precision,      intent(in)    :: vor_bar(3)
            integer                              :: n
            double precision                     :: dvor(3), asum
            double precision                     :: buf(4)

            asum = zero
            dvor = zero

            !$omp parallel default(shared)
            !$omp do private(n) reduction(+: dvor, asum)
            do n = 1, pcont%local_num
                asum = asum + pcont%area(n)
                dvor = dvor + pcont%vorticity(:, n) * pcont%area(n)
            enddo
            !$omp end do
            !$omp end parallel

            buf(1) = asum
            buf(2:4) = dvor
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:4),                &
                               4,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of parcel_correction::m_apply_surface_vortcor.")
            asum = buf(1)
            dvor = buf(2:4)

            ! vorticity correction:
            ! vor_mean = dvor / asum
            ! if vor_mean < vor_bar --> correction < 0 --> we must add to parcel vorticity
            ! if vor_mean > vor_bar --> correction > 0 --> we must subtract from parcel vorticity
            dvor = dvor / asum - vor_bar

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, pcont%local_num
                pcont%vorticity(:, n) = pcont%vorticity(:, n) - dvor
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine m_apply_surface_vortcor

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Ensure that the total parcel area sum(A_i) equals the domain area
    ! (Lx * Ly = extent(1)*extent(2))
    subroutine area_correction
        double precision :: total_area(2)

        total_area(1) = sum(bot_parcels%area(1:bot_parcels%local_num))
        total_area(2) = sum(top_parcels%area(1:top_parcels%local_num))

        call MPI_Allreduce(MPI_IN_PLACE,            &
                           total_area(1:2),         &
                           2,                       &
                           MPI_DOUBLE_PRECISION,    &
                           MPI_SUM,                 &
                           world%comm,              &
                           world%err)

        call mpi_check_for_error(world, &
            "in MPI_Allreduce of parcel_correction::area_correction.")

        ! Correction factor: f = Lx Ly / sum Ai.
        total_area = extent(1) * extent(2) / total_area

        bot_parcels%area(1:bot_parcels%local_num) = bot_parcels%area(1:bot_parcels%local_num) * total_area(1)
        top_parcels%area(1:top_parcels%local_num) = top_parcels%area(1:top_parcels%local_num) * total_area(2)

    end subroutine area_correction

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_surf_laplace

        call area2grid(l_halo_swap=.true.)

        call m_apply_laplace(0,  bot_parcels)
        call m_apply_laplace(nz, top_parcels)
    end subroutine apply_surf_laplace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_apply_laplace(iz, surf_parcels)
        integer,               intent(in)    :: iz
        type(ellipse_pc_type), intent(inout) :: surf_parcels
        double precision                     :: phi( box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision                     :: grad(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 2)
        double precision                     :: weights(0:1, 0:1)
        integer                              :: n, l, is, js

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(iz, :, :) * acelli - one

        call m_diverge(phi, grad)

        !------------------------------------------------------------------
        ! Increment parcel positions using (ud, vd) field:
        !$omp parallel default(shared)
        !$omp do private(n, l, is, js, weights)
        do n = 1, surf_parcels%local_num
            call bilinear(surf_parcels%position(:, n), is, js, weights)

            do l = 1, 2
                surf_parcels%position(l, n) = surf_parcels%position(l, n)               &
                                            + sum(weights * grad(js:js+1, is:is+1, l))
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        call parcel_communicate(surf_parcels)

    end subroutine m_apply_laplace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_surf_gradient(prefactor, max_compression)
        double precision, intent(in) :: prefactor
        double precision, intent(in) :: max_compression

        call area2grid(l_halo_swap=.true.)

        call m_apply_gradient(0,  bot_parcels, prefactor, max_compression)
        call m_apply_gradient(nz, top_parcels, prefactor, max_compression)

    end subroutine apply_surf_gradient

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_apply_gradient(iz, surf_parcels, prefactor, max_compression)
        integer,               intent(in)    :: iz
        type(ellipse_pc_type), intent(inout) :: surf_parcels
        double precision,      intent(in)    :: prefactor
        double precision,      intent(in)    :: max_compression
        double precision                     :: phi(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision                     :: weights(0:1, 0:1)
        double precision                     :: xs, ys, xf, yf, xfc, yfc, lim_x, lim_y
        integer                              :: n, is, js

        ! form divergence field * dt and store in phi temporarily:
        phi = volg(iz, :, :) * acelli - one

        !$omp parallel default(shared)
        !$omp do private(n, is, js, weights, xs, ys, xf, yf, xfc, yfc, lim_x, lim_y)
        do n = 1, surf_parcels%local_num

            call bilinear(surf_parcels%position(:, n), is, js, weights)

            xf = weights(0, 1) + weights(1, 1) ! fractional position along x
            yf = weights(1, 0) + weights(1, 1) ! fractional position along y
            xfc = one - xf
            yfc = one - yf


            lim_x = dx(1) * xf * xfc
            xs = - prefactor * lim_x * (                             &
                            yfc * (phi(js,   is+1) - phi(js,   is))  &
                          + yf  * (phi(js+1, is+1) - phi(js+1, is)))

            lim_x = lim_x * max_compression
            xs = max(-lim_x, min(xs, lim_x))

            lim_y = dx(2) * yf * yfc
            ys = - prefactor * lim_y * (                            &
                            xfc * (phi(js+1, is  ) - phi(js, is))   &
                          + xf  * (phi(js+1, is+1) - phi(js, is+1)))

            lim_y = lim_y * max_compression
            ys = max(-lim_y, min(ys, lim_y))

            surf_parcels%position(1, n) = surf_parcels%position(1, n) + xs
            surf_parcels%position(2, n) = surf_parcels%position(2, n) + ys
        enddo
        !$omp end do
        !$omp end parallel

        call parcel_communicate(surf_parcels)

    end subroutine m_apply_gradient

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Computes a divergent flow field (ud, vd) = grad(phi) where
    ! Lap(phi) = div (given).
    subroutine m_diverge(div, grad)
        double precision, intent(inout) :: div(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
        double precision, intent(out)   :: grad(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), 2)
        double precision                :: ds(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                :: us(box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                                           vs(box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        !------------------------------------------------------------------
        ! Convert phi to spectral space as ds:
        call fft2d(div, ds)

        ! Invert Laplace's operator spectrally
        call lapinv(ds)

        ! Compute x derivative spectrally:
        call diff2dx(ds, us)

        ! Reverse FFT to define x velocity component ud:
        call ifft2d(us, grad(:, :, 1))

        ! Compute y derivative spectrally:
        call diff2dy(ds, vs)

        ! Reverse FFT to define y velocity component vd:
        call ifft2d(vs, grad(:, :, 2))

        ! Fill halo grid points:
        call field_halo_fill_vector(grad, l_alloc=.true.)

    end subroutine m_diverge

end module surface_parcel_correction
