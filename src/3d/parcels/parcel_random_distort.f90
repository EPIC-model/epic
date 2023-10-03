module parcel_random_distort
    use constants, only : pi, two
    use mpi_timer, only : start_timer, stop_timer
    use parameters, only : dx
    use parcel_container, only : parcels, n_parcels
    use fields
    use omp_lib
    use mpi_layout
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use parcel_bc, only : apply_parcel_reflective_bc
    use parcel_mpi, only : parcel_communicate
    use parcel_interpl, only : trilinear
    use parcel_ellipsoid, only : get_B33, I_B11, I_B12, I_B13, I_B22, I_B23
    implicit none

    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)
    double precision :: random_strain(9)
    integer, parameter :: I_RAND_DUDX=1
    integer, parameter :: I_RAND_DUDY=2
    integer, parameter :: I_RAND_DUDZ=3
    integer, parameter :: I_RAND_DVDX=4
    integer, parameter :: I_RAND_DVDY=5
    integer, parameter :: I_RAND_DVDZ=6
    integer, parameter :: I_RAND_DWDX=7
    integer, parameter :: I_RAND_DWDY=8
    integer, parameter :: I_RAND_DWDZ=9


    public :: parcel_stochastic_distort

    contains

        subroutine parcel_stochastic_distort(dt)
            double precision, intent(in)  :: dt
            double precision  :: prefac, u_random, v_random, p_strain, x
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights, p_strain, u_random, v_random, random_strain, prefac, x)
            do n = 1, n_parcels
                call trilinear(parcels%position(:, n), is, js, ks, weights)
                p_strain=sum(weights * strain_mag(ks:ks+1, js:js+1, is:is+1))
                prefac=0.5*sqrt(p_strain*dt)*p_strain
                call gauss_bound_random_number(x,3.0d0)
                u_random=sqrt(1.5)*x
                call gauss_bound_random_number(x,3.0d0)
                v_random=sqrt(1.5)*x

                ! First, sort out diagonal
                random_strain(I_RAND_DUDX)=prefac*((1./sqrt(2.))*u_random+(1./sqrt(6.))*v_random)
                random_strain(I_RAND_DVDY)=prefac*(-(1./sqrt(2.))*u_random+(1./sqrt(6.))*v_random)
                random_strain(I_RAND_DWDZ)=prefac*(-random_strain(I_RAND_DUDX)-random_strain(I_RAND_DVDY))

                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DUDY)=prefac*x
                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DUDZ)=prefac*x
                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DVDX)=prefac*x
                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DVDZ)=prefac*x
                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DWDX)=prefac*x
                call gauss_bound_random_number(x,3.0d0)
                random_strain(I_RAND_DWDY)=prefac*x

                parcels%B(:, n) = parcels%B(:, n)               &
                                      + dt*get_dBdt_random_matrix(parcels%B(:, n),           &
                                                 random_strain(:),      &
                                                 parcels%volume(n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine parcel_stochastic_distort

       ! Advance the B matrix.
        ! @param[in] Bin are the B matrix components of the parcel
        ! @param[in] S is the local velocity strain
        ! @param[in] vorticity of parcel
        ! @param[in] volume is the parcel volume
        ! @returns dB/dt in Bout
        function get_dBdt_random_matrix(Bin, random_strain, volume) result(Bout)
            double precision, intent(in) :: Bin(I_B23)
            double precision, intent(in) :: random_strain(9)
            double precision, intent(in) :: volume
            double precision             :: Bout(5), B33

            B33 = get_B33(Bin, volume)

            ! dB11/dt = 2 * (du/dx * B11 + du/dy * B12 + du/dz * B13)
            Bout(I_B11) = two * (random_strain(I_RAND_DUDX) * Bin(I_B11) + &
                          random_strain(I_RAND_DUDY) * Bin(I_B12) + &
                          random_strain(I_RAND_DUDZ) * Bin(I_B13))

            ! dB12/dt =
            Bout(I_B12) = random_strain(I_RAND_DVDX) * Bin(I_B11) & !   dv/dx * B11
                        - random_strain(I_RAND_DWDZ) * Bin(I_B12) & ! - dw/dz * B12
                        + random_strain(I_RAND_DVDZ) * Bin(I_B13) & ! + dv/dz * B13
                        + random_strain(I_RAND_DUDY) * Bin(I_B22) & ! + du/dy * B22
                        + random_strain(I_RAND_DUDZ) * Bin(I_B23)   ! + du/dz * B23

            ! dB13/dt =
            Bout(I_B13) = random_strain(I_RAND_DWDX) * Bin(I_B11) & !   dw/dx * B11
                        + random_strain(I_RAND_DWDY) * Bin(I_B12) & ! + dw/dy * B12
                        - random_strain(I_RAND_DVDY) * Bin(I_B13) & ! - dv/dy * B13
                        + random_strain(I_RAND_DUDY) * Bin(I_B23) & ! + du/dy * B23
                        + random_strain(I_RAND_DUDZ) * B33          ! + du/dz * B33

            ! dB22/dt = 2 * (dv/dx * B12 + dv/dy * B22 + dv/dz * B23)
            Bout(I_B22) = two * (random_strain(I_RAND_DVDX) * Bin(I_B12) + &
                          random_strain(I_RAND_DVDY) * Bin(I_B22) +  &
                          random_strain(I_RAND_DVDZ) * Bin(I_B23))

            ! dB23/dt =
            Bout(I_B23) =  random_strain(I_RAND_DWDX) * Bin(I_B12) & !   dw/dx * B12
                        +  random_strain(I_RAND_DVDX) * Bin(I_B13) & ! + dv/dx * B13
                        +  random_strain(I_RAND_DWDY) * Bin(I_B22) & ! + dw/dy * B22
                        -  random_strain(I_RAND_DUDX) * Bin(I_B23) & ! - du/dx * B23
                        +  random_strain(I_RAND_DVDZ) * B33          ! + dv/dz * B33
        end function get_dBdt_random_matrix

        subroutine gauss_bound_random_number(x,bound)
            implicit none
            double precision,intent(out) :: x
            double precision,intent(in) :: bound
            double precision :: r1,r2
            call random_number(r1)
            call random_number(r2)
            x = sqrt(-2.0*log(r1))*cos(2.0*pi*r2)
            x = max(min(x,bound),-bound)
        end subroutine gauss_bound_random_number

end module parcel_random_distort
