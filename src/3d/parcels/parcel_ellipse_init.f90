submodule (parcel_ellipse) parcel_ellipse_init_smod
    use options, only : parcel
    use constants, only : one
    use parameters, only : nx, ny, acell, dx, dxi, lower, vcell, max_num_surf_parcels
    use mpi_environment, only : world, MPI_SUM
    use mpi_layout, only : box
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    contains

        module subroutine parcel_ellipse_init(this)
            class(ellipse_pc_type), intent(inout) :: this
            integer                               :: ix, i, iy, j, k, n_per_dim, n
            double precision                      :: im, corner(2), a2, lam, ratio

            call this%allocate(max_num_surf_parcels)

            !------------------------------------------------------------------
            ! Set the number of parcels
            this%local_num = parcel%n_surf_per_cell * box%size(1) * box%size(2)

            if (this%local_num > this%max_num) then
                print *, "Number of parcels exceeds limit of", &
                          this%max_num, ". Exiting."
                call mpi_exit_on_error
            endif

            !------------------------------------------------------------------
            ! Initialise parcel positions on a regular grid:

            this%total_num = this%local_num
            if (world%size > 1) then
                call mpi_blocking_reduce(this%total_num, MPI_SUM, world)
            endif

            ! number of parcels per dimension
            n_per_dim = int(dsqrt(dble(parcel%n_surf_per_cell)))
            if (n_per_dim ** 2 .ne. parcel%n_surf_per_cell) then
                if (world%rank == world%root) then
                    print *, "Number of parcels per cell (", &
                             parcel%n_surf_per_cell, ") not a square."
                endif
                call mpi_exit_on_error
            endif

            im = one / dble(n_per_dim)

            k = 1
            do iy = box%lo(2), box%hi(2)
                do ix = box%lo(1), box%hi(1)
                    corner = lower(1:2) + dble((/ix, iy/)) * dx(1:2)
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            this%position(1, k) = corner(1) + dx(1) * (dble(i) - f12) * im
                            this%position(2, k) = corner(2) + dx(2) * (dble(j) - f12) * im
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo

            if (.not. this%local_num == k - 1) then
                call mpi_exit_on_error(&
                    "in parcel_ellipse_init: Number of parcels disagree!")
            endif

            !------------------------------------------------------------------
            ! Initialize the area of each parcel:
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, this%local_num
                this%area(n) = acell / dble(parcel%n_surf_per_cell)
            enddo
            !$omp end do
            !$omp end parallel


            !------------------------------------------------------------------
            ! Initialise the shape matrix of each parcel:
            ratio = dx(1) / dx(2)

            ! aspect ratio: lam = a / b
            lam = max(dx(2) / dx(1), ratio)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, this%local_num
                ! B11
                this%B(1, n) = ratio * this%get_ab(this%area(n))

                ! B12
                this%B(2, n) = zero

                ! B22
                this%B(3, n) = this%get_ab(this%area(n)) / ratio
            enddo
            !$omp end do
            !$omp end parallel

            !------------------------------------------------------------------
            ! Ensure the initial parcels are not extremely elongated:
            ! --> refine by splitting
            do while (lam >= parcel%lambda_max)
                call this%split
                ! all parcels have the same initials aspect ratio, we
                ! need to check only one of them
                a2 = this%get_eigenvalue(1)
                lam = a2 / this%get_ab(this%area(1))
            end do

            !------------------------------------------------------------------
            ! Initialise all other parcel attributes to their default values:
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, this%local_num
                this%vorticity(:, n) = zero
                this%buoyancy(n)     = zero
#ifndef ENABLE_DRY_MODE
                this%humidity(n)     = zero
#endif
                this%volume(n) = vcell / dble(parcel%n_surf_per_cell)
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine parcel_ellipse_init

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end submodule parcel_ellipse_init_smod
