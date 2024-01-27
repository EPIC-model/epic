! =============================================================================
!     Parcel ellipsoid container for 3D parcels.
!
!     Reference:
!       Dritschel, D., Reinaud, J., & McKiver, W. (2004).
!       The quasi-geostrophic ellipsoidal vortex model.
!       Journal of Fluid Mechanics, 505, 201-223.
!       doi:10.1017/S0022112004008377
! =============================================================================
module parcel_ellipsoid
    use parcel_container, only : pc_type
    use parameters, only : nz
    use dimensions, only : I_X, I_Y, I_Z
    use constants, only : fpi   &
                        , fpi4  &
                        , f34   &
                        , zero  &
                        , two   &
                        , three &
                        , five  &
                        , seven
    use scherzinger, only : scherzinger_diagonalise &
                          , scherzinger_eigenvalues
    use mpi_utils, only : mpi_exit_on_error
    use armanip, only : resize_array
    implicit none

    double precision, parameter :: rho = dsqrt(two / five)
    double precision, parameter :: f3pi4 = three * fpi4
    double precision, parameter :: f5pi4 = five * fpi4
    double precision, parameter :: f7pi4 = seven * fpi4
    double precision, parameter :: costheta(4) = dcos((/fpi4, f3pi4, f5pi4, f7pi4/))
    double precision, parameter :: sintheta(4) = dsin((/fpi4, f3pi4, f5pi4, f7pi4/))

    type, extends(pc_type) :: ellipsoid_pc_type

        ! ---------------------------------------------------------------------
        ! Additional parcel attributes:
        double precision, allocatable :: Vetas(:, :), Vtaus(:, :)

        integer :: IDX_ELL_VETA, IDX_ELL_VTAU

        contains

            procedure          :: setup => parcel_ellipsoid_setup
            procedure          :: dealloc => parcel_ellipsoid_dealloc
            procedure          :: replace => parcel_ellipsoid_replace
            procedure          :: resize => parcel_ellipsoid_resize
            procedure          :: get_points => parcel_ellipsoid_get_points
            procedure, private :: get_full_matrix => parcel_ellipsoid_get_full_matrix
            procedure          :: get_eigenvalues => parcel_ellipsoid_get_eigenvalues
            procedure          :: get_determinant => parcel_ellipsoid_get_determinant
            procedure, private :: get_eigenvectors => parcel_ellipsoid_get_eigenvectors
            procedure          :: diagonalise => parcel_ellipsoid_diagonalise
            procedure          :: get_B33 => parcel_ellipsoid_get_B33
            procedure          :: get_abc => parcel_ellipsoid_get_abc
            procedure          :: get_aspect_ratio => parcel_ellipsoid_get_aspect_ratio
            procedure, private :: get_angles => parcel_ellipsoid_get_angles
            procedure          :: get_local_cell_index => parcel_ellipsoid_get_local_cell_index
            procedure          :: is_small => parcel_ellipsoid_is_small

    end type ellipsoid_pc_type

    !--------------------------------------------------------------------------
    ! Define interface for submodule routines:
    interface
        ! Implemented in parcel_ellipse_init
        module subroutine parcel_ellipsoid_get_local_cell_index(this, n, ix, iy, iz)
            class(ellipsoid_pc_type), intent(in)  :: this
            integer,                  intent(in)  :: n
            integer,                  intent(out) :: ix, iy, iz
        end subroutine parcel_ellipsoid_get_local_cell_index

        module logical function parcel_ellipsoid_is_small(this, n)
            class(ellipsoid_pc_type), intent(in)  :: this
            integer,                  intent(in)  :: n
        end function parcel_ellipsoid_is_small

    end  interface


    integer, parameter :: I_B11 = 1 & ! index for B11 matrix component
                        , I_B12 = 2 & ! index for B12 matrix component
                        , I_B13 = 3 & ! index for B13 matrix component
                        , I_B22 = 4 & ! index for B22 matrix component
                        , I_B23 = 5   ! index for B23 matrix component


    private :: rho                  &
             , f3pi4                &
             , f5pi4                &
             , f7pi4                &
             , costheta             &
             , sintheta

    contains

        ! Sets the buffer indices for sending parcels around
        subroutine parcel_ellipsoid_setup(this, num)
            class(ellipsoid_pc_type), intent(inout) :: this
            integer,                  intent(in)    :: num
            integer                                 :: i

            call this%alloc(num=num,    &
                            n_pos=3,    &
                            n_vec=3,    &
                            n_shape=5,  &
                            n_strain=5)

            allocate(this%Vetas(3, num))
            allocate(this%Vtaus(3, num))

            this%IDX_POS_BEG   = 1   ! x-position
            this%IDX_POS_END   = 3   ! z-position
            this%IDX_VOR_BEG   = 4   ! x-vorticity
            this%IDX_VOR_END   = 6   ! z-vorticity
            this%IDX_SHAPE_BEG = 7   ! B11 shape matrix element
            this%IDX_SHAPE_END = 11  ! B23 shape matrix element
            this%IDX_VOL       = 12  ! volume
            this%IDX_BUO       = 13  ! buoyancy

            i = this%IDX_BUO + 1
#ifndef ENABLE_DRY_MODE
            this%IDX_HUM  = i
            i = i + 1
#endif

            ! LS-RK variables
            this%IDX_RK_POS_BEG    = i
            this%IDX_RK_POS_END    = i + 2
            this%IDX_RK_VOR_BEG    = i + 3
            this%IDX_RK_VOR_END    = i + 5
            this%IDX_RK_SHAPE_BEG  = i + 6
            this%IDX_RK_SHAPE_END  = i + 10
            this%IDX_RK_STRAIN_BEG = i + 11
            this%IDX_RK_STRAIN_END = i + 15

            this%IDX_ELL_VETA = i + 16
            this%IDX_ELL_VTAU = i + 19

            this%attr_num = this%IDX_ELL_VTAU + 2

            ! these parcels are in all vertical grid cells
            this%nz = nz

        end subroutine parcel_ellipsoid_setup

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipsoid_dealloc(this)
            class(ellipsoid_pc_type), intent(inout) :: this

            call this%dealloc

            if (.not. allocated(this%Vetas)) then
                return
            endif

            deallocate(this%Vetas)
            deallocate(this%Vtaus)

        end subroutine parcel_ellipsoid_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipsoid_replace(this, n, m)
            class(ellipsoid_pc_type), intent(inout) :: this
            integer,                  intent(in)    :: n, m

            ! Call parent class subroutine
            call this%replace(n, m)

            this%Vetas(:, n) = this%Vetas(:, m)
            this%Vtaus(:, n) = this%Vtaus(:, m)

        end subroutine parcel_ellipsoid_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipsoid_resize(this, new_size)
            class(ellipsoid_pc_type), intent(inout) :: this
            integer,                  intent(in)    :: new_size

            ! Call parent class subroutine
            call this%resize(new_size)

            call resize_array(this%Vetas, new_size, n_copy=this%local_num)
            call resize_array(this%Vtaus, new_size, n_copy=this%local_num)

        end subroutine parcel_ellipsoid_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipsoid_serialize(this, n, buffer)
            class(ellipsoid_pc_type), intent(in)  :: this
            integer,                  intent(in)  :: n
            double precision,         intent(out) :: buffer(this%attr_num)

            ! Call parent class subroutine
            call this%serialize(n, buffer)

            buffer(this%IDX_ELL_VETA:this%IDX_ELL_VETA+2) = this%Vetas(:, n)
            buffer(this%IDX_ELL_VTAU:this%IDX_ELL_VTAU+2) = this%Vtaus(:, n)

        end subroutine parcel_ellipsoid_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipsoid_deserialize(this, n, buffer)
            class(ellipsoid_pc_type), intent(inout) :: this
            integer,                  intent(in)    :: n
            double precision,         intent(in)    :: buffer(this%attr_num)

            ! Call parent class subroutine
            call this%deserialize(n, buffer)

            this%Vetas(:, n) = buffer(this%IDX_ELL_VETA:this%IDX_ELL_VETA+2)
            this%Vtaus(:, n) = buffer(this%IDX_ELL_VTAU:this%IDX_ELL_VTAU+2)

        end subroutine parcel_ellipsoid_deserialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the ellipsoid support points for par2grid and grid2par
        ! @param[in] n parcel index
        ! @param[in] l_reuse (optional) if support points should be reused and not calculated.
        ! @returns the parcel support points
        function parcel_ellipsoid_get_points(this, n, l_reuse) result(points)
            class(ellipsoid_pc_type), intent(inout) :: this
            integer,                  intent(in)    :: n
            logical, optional,        intent(in)    :: l_reuse
            double precision                        :: Veta(3), Vtau(3), D(3), V(3, 3)
            integer                                 :: j
            double precision                        :: points(3, 4)


            if (present(l_reuse)) then
                if (l_reuse) then
                    Veta = this%Vetas(:, n)
                    Vtau = this%Vtaus(:, n)
                else
                    call this%diagonalise(n, D, V)
                    Veta = dsqrt(dabs(D(I_X) - D(I_Z))) * rho * V(:, I_X)
                    Vtau = dsqrt(dabs(D(I_Y) - D(I_Z))) * rho * V(:, I_Y)

                    this%Vetas(:, n) = Veta
                    this%Vtaus(:, n) = Vtau
                endif
            else
                ! (/a2, b2, c2/) with a >= b >= c
                ! D = (/a2, b2, c2/)
                call this%diagonalise(n, D, V)

                Veta = dsqrt(dabs(D(I_X) - D(I_Z))) * rho * V(:, I_X)
                Vtau = dsqrt(dabs(D(I_Y) - D(I_Z))) * rho * V(:, I_Y)

                this%Vetas(:, n) = Veta
                this%Vtaus(:, n) = Vtau
            endif

            do j = 1, 4
                ! support point in the reference frame of the ellipsoid
                ! theta = j * pi / 2 - pi / 4 (j = 1, 2, 3, 4)
                ! x_j = eta * rho * cos(theta_j)
                ! y_j = tau * rho * sin(theta_j)

                ! suppport point in the global reference frame
                points(:, j) = this%position(:, n) &
                             + Veta * costheta(j)  &
                             + Vtau * sintheta(j)
            enddo
        end function parcel_ellipsoid_get_points

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the parcel shape matrix.
        ! @param[in] n parcel index
        ! @returns the upper trinagular matrix
        function parcel_ellipsoid_get_full_matrix(this, n) result(U)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: U(3, 3)

            U(1, 1) = this%B(I_B11, n)
            U(1, 2) = this%B(I_B12, n)
            U(1, 3) = this%B(I_B13, n)
            U(2, 1) = U(1, 2)
            U(2, 2) = this%B(I_B22, n)
            U(2, 3) = this%B(I_B23, n)
            U(3, 1) = U(1, 3)
            U(3, 2) = U(2, 3)
            U(3, 3) = this%get_B33(n)
        end function parcel_ellipsoid_get_full_matrix

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain all eigenvalues sorted in descending order
        ! @param[in] n parcel index
        ! @returns all eigenvalues (sorted in descending order)
        function parcel_ellipsoid_get_eigenvalues(this, n) result(D)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: U(3, 3)
            double precision                     :: D(3)

            U = this%get_full_matrix(n)

            call scherzinger_eigenvalues(U, D)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                call mpi_exit_on_error(&
                    "in parcel_ellipsoid::get_eigenvalues: Invalid parcel shape.")
            endif
#endif
        end function parcel_ellipsoid_get_eigenvalues

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the determinant of the shape matrix
        ! @param[in] n parcel index
        ! @returns the determinant det(B)
        function parcel_ellipsoid_get_determinant(this, n) result(detB)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: detB, B33

            B33 = this%get_B33(n)

            detB = this%B(I_B11, n) * (this%B(I_B22, n) * B33              - this%B(I_B23, n) ** 2)               &
                 - this%B(I_B12, n) * (this%B(I_B12, n) * B33              - this%B(I_B13, n) * this%B(I_B23, n)) &
                 + this%B(I_B13, n) * (this%B(I_B12, n) * this%B(I_B23, n) - this%B(I_B13, n) * this%B(I_B22, n))
        end function parcel_ellipsoid_get_determinant

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the normalized eigenvectors.
        ! The eigenvector V(:, j) belongs to the j-th
        ! eigenvalue.
        ! @param[in] n parcel index
        ! @returns the eigenvectors
        function parcel_ellipsoid_get_eigenvectors(this, n) result(V)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: U(3, 3), D(3), V(3, 3)

            U = this%get_full_matrix(n)

            call scherzinger_diagonalise(U, D, V)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                call mpi_exit_on_error(&
                    "in parcel_ellipsoid::get_eigenvectors: Invalid parcel shape.")
            endif
#endif
        end function parcel_ellipsoid_get_eigenvectors

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the eigenvalue decomposition B = V^T * D * V
        ! where D has the eigenvalues on its diagonal
        ! and V contains the eigenvectors in its columns.
        ! The eigenvector V(:, j) belongs to the j-th
        ! eigenvalue.
        ! @param[in] n parcel index
        ! @param[out] D eigenvalues (sorted in descending order)
        ! @param[out] V eigenvectors
        subroutine parcel_ellipsoid_diagonalise(this, n, D, V)
            class(ellipsoid_pc_type), intent(in)  :: this
            integer,                  intent(in)  :: n
            double precision,         intent(out) :: D(3), V(3, 3)
            double precision                      :: U(3, 3)

            U = this%get_full_matrix(n)

            call scherzinger_diagonalise(U, D, V)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                call mpi_exit_on_error(&
                    "in parcel_ellipsoid::diagonalise: Invalid parcel shape.")
            endif
#endif
        end subroutine parcel_ellipsoid_diagonalise

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the B33 matrix element
        ! @param[in] n parcel index
        ! @returns B33
        function parcel_ellipsoid_get_B33(this, n) result(B33)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: abc
            double precision                     :: B33

            abc = this%get_abc(this%volume(n))

            if (dabs(this%B(I_B11, n) * this%B(I_B22, n) - this%B(I_B12, n) ** 2) <= epsilon(abc)) then
                call mpi_exit_on_error(&
                    "in parcel_ellipsoid::get_B33: Division by small number!")
            endif

            B33 = (abc ** 2 - this%B(I_B13, n) * (this%B(I_B12, n) * this%B(I_B23, n) -     &
                                                  this%B(I_B13, n) * this%B(I_B22, n))      &
                            + this%B(I_B11, n) * this%B(I_B23, n) ** 2                      &
                            - this%B(I_B12, n) * this%B(I_B13, n) * this%B(I_B23, n))       &
                / (this%B(I_B11, n) * this%B(I_B22, n) - this%B(I_B12, n) ** 2)

        end function parcel_ellipsoid_get_B33

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] volume
        ! @returns abc = 3 * volume / (4 * pi)
        pure elemental function parcel_ellipsoid_get_abc(this, volume) result(abc)
            class(ellipsoid_pc_type), intent(in) :: this
            double precision,         intent(in) :: volume
            double precision                     :: abc

            abc = f34 * volume * fpi
        end function parcel_ellipsoid_get_abc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the aspect ratio a/c of the parcel(s).
        ! @param[in] D eigenvalues sorted in descending order
        ! @param[in] volume of the parcel(s)
        ! @returns a/c
        pure function parcel_ellipsoid_get_aspect_ratio(this, D) result(lam)
            class(ellipsoid_pc_type), intent(in) :: this
            double precision,         intent(in) :: D(3)
            double precision                     :: lam

            lam = dsqrt(D(I_X) / D(I_Z))
        end function parcel_ellipsoid_get_aspect_ratio

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! @param[in] n parcel index
        function parcel_ellipsoid_get_angles(this, n) result(angles)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            double precision                     :: evec(3, 3)
            double precision                     :: angles(2) ! (/azimuth, polar/)

            evec = this%get_eigenvectors(n)

            ! azimuthal angle
            angles(I_X) = datan2(evec(I_Y, I_X), evec(I_X, I_X))

            ! polar angle
            angles(I_Y) = dasin(evec(I_Z, I_Z))

        end function parcel_ellipsoid_get_angles

end module parcel_ellipsoid
