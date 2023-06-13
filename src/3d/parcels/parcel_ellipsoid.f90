! =============================================================================
!     This module is the 3D version and contains all ellipsoid operations.
!     Reference:
!       Dritschel, D., Reinaud, J., & McKiver, W. (2004).
!       The quasi-geostrophic ellipsoidal vortex model.
!       Journal of Fluid Mechanics, 505, 201-223.
!       doi:10.1017/S0022112004008377
! =============================================================================
module parcel_ellipsoid
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use constants, only : fpi   &
                        , fpi4  &
                        , f34   &
                        , zero  &
                        , two   &
                        , three &
                        , five  &
                        , seven
    use parameters, only : max_num_parcels
    use scherzinger, only : scherzinger_diagonalise &
                          , scherzinger_eigenvalues
    implicit none

    double precision, parameter :: rho = dsqrt(two / five)
    double precision, parameter :: f3pi4 = three * fpi4
    double precision, parameter :: f5pi4 = five * fpi4
    double precision, parameter :: f7pi4 = seven * fpi4
    double precision, parameter :: costheta(4) = dcos((/fpi4, f3pi4, f5pi4, f7pi4/))
    double precision, parameter :: sintheta(4) = dsin((/fpi4, f3pi4, f5pi4, f7pi4/))

    double precision, allocatable :: Vetas(:, :), Vtaus(:, :)

    integer, parameter :: I_B11 = 1 & ! index for B11 matrix component
                        , I_B12 = 2 & ! index for B12 matrix component
                        , I_B13 = 3 & ! index for B13 matrix component
                        , I_B22 = 4 & ! index for B22 matrix component
                        , I_B23 = 5   ! index for B23 matrix component

    private :: rho, f3pi4, f5pi4, f7pi4, costheta, sintheta, get_upper_triangular, Vetas, Vtaus

    contains

        subroutine parcel_ellipsoid_allocate
            allocate(Vetas(n_dim, max_num_parcels))
            allocate(Vtaus(n_dim, max_num_parcels))
        end subroutine parcel_ellipsoid_allocate

        subroutine parcel_ellipsoid_deallocate
            deallocate(Vetas)
            deallocate(Vtaus)
        end subroutine parcel_ellipsoid_deallocate

        ! Obtain the parcel shape matrix.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns the upper trinagular matrix
        function get_upper_triangular(B, volume) result(U)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(n_dim, n_dim)

            U(1, 1) = B(I_B11)
            U(1, 2) = B(I_B12)
            U(1, 3) = B(I_B13)
            U(2, 2) = B(I_B22)
            U(2, 3) = B(I_B23)
            U(3, 3) = get_B33(B, volume)
        end function get_upper_triangular

        ! Obtain all eigenvalues sorted in descending order
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns all eigenvalues (sorted in descending order)
        function get_eigenvalues(B, volume) result(D)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(n_dim, n_dim)
            double precision             :: D(n_dim)

            U = get_upper_triangular(B, volume)

            call scherzinger_eigenvalues(U, D)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                print *, "Invalid parcel shape."
                stop
            endif
#endif
        end function get_eigenvalues

        function get_determinant(B, volume) result(detB)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: detB, B33

            B33 = get_B33(B, volume)

            detB = B(I_B11) * (B(I_B22) * B33 - B(I_B23) ** 2)            &
                 - B(I_B12) * (B(I_B12) * B33 - B(I_B13) * B(I_B23))      &
                 + B(I_B13) * (B(I_B12) * B(I_B23) - B(I_B13) * B(I_B22))
        end function get_determinant

        ! Obtain the normalized eigenvectors.
        ! The eigenvector V(:, j) belongs to the j-th
        ! eigenvalue.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns the eigenvectors
        function get_eigenvectors(B, volume) result(V)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(n_dim, n_dim), D(n_dim), V(n_dim, n_dim)

            U = get_upper_triangular(B, volume)

            call scherzinger_diagonalise(U, D, V)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                print *, "Invalid parcel shape."
                stop
            endif
#endif
        end function get_eigenvectors

        ! Compute the eigenvalue decomposition B = V^T * D * V
        ! where D has the eigenvalues on its diagonal
        ! and V contains the eigenvectors in its columns.
        ! The eigenvector V(:, j) belongs to the j-th
        ! eigenvalue.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @param[out] D eigenvalues (sorted in descending order)
        ! @param[out] V eigenvectors
        subroutine diagonalise(B, volume, D, V)
            double precision, intent(in)  :: B(5)
            double precision, intent(in)  :: volume
            double precision, intent(out) :: D(n_dim), V(n_dim, n_dim)
            double precision              :: U(n_dim, n_dim)

            U = get_upper_triangular(B, volume)

            call scherzinger_diagonalise(U, D, V)

#ifndef NDEBUG
            ! check if any eigenvalue is less or equal zero
            if (minval(D) <= zero) then
                print *, "Invalid parcel shape."
                stop
            endif
#endif
        end subroutine diagonalise

        ! Obtain the B33 matrix element
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns B33
        function get_B33(B, volume) result(B33)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: abc
            double precision             :: B33

            abc = get_abc(volume)

            if (dabs(B(I_B11) * B(I_B22) - B(I_B12) ** 2) <= epsilon(abc)) then
                print *, "Error in get_B33: Division by small number!"
                stop
            endif

            B33 = (abc ** 2 - B(I_B13) * (B(I_B12) * B(I_B23) - B(I_B13) * B(I_B22)) &
                            + B(I_B11) * B(I_B23) ** 2                               &
                            - B(I_B12) * B(I_B13) * B(I_B23))                        &
                / (B(I_B11) * B(I_B22) - B(I_B12) ** 2)

        end function get_B33

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] volume of the parcel
        ! @returns abc = 3 * volume / (4 * pi)
        pure elemental function get_abc(volume) result(abc)
            double precision, intent(in) :: volume
            double precision             :: abc

            abc = f34 * volume * fpi
        end function get_abc

        ! Obtain the aspect ratio a/c of the parcel(s).
        ! @param[in] D eigenvalues sorted in descending order
        ! @param[in] volume of the parcel(s)
        ! @returns a/c
        pure function get_aspect_ratio(D) result(lam)
            double precision, intent(in) :: D(n_dim)
            double precision             :: lam

            lam = dsqrt(D(I_X) / D(I_Z))
        end function get_aspect_ratio

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] position vector of the parcel
        ! @param[in] volume of the parcel
        ! @param[in] B matrix elements of the parcel
        ! @returns the parcel support points
        function get_ellipsoid_points(position, volume, B, n, l_reuse) result(points)
            double precision,  intent(in) :: position(n_dim)
            double precision,  intent(in) :: volume
            double precision,  intent(in) :: B(5)        ! B11, B12, B13, B22, B23
            integer, optional, intent(in) :: n
            logical, optional, intent(in) :: l_reuse
            double precision              :: Veta(n_dim), Vtau(n_dim), D(n_dim), V(n_dim, n_dim)
            integer                       :: j
            double precision              :: points(n_dim, 4)


            if (present(l_reuse)) then
                if (l_reuse) then
                    Veta = Vetas(:, n)
                    Vtau = Vtaus(:, n)
                else
                    call diagonalise(B, volume, D, V)
                    Veta = dsqrt(dabs(D(I_X) - D(I_Z))) * rho * V(:, I_X)
                    Vtau = dsqrt(dabs(D(I_Y) - D(I_Z))) * rho * V(:, I_Y)

                    Vetas(:, n) = Veta
                    Vtaus(:, n) = Vtau
                endif
            else
                ! (/a2, b2, c2/) with a >= b >= c
                ! D = (/a2, b2, c2/)
                call diagonalise(B, volume, D, V)

                Veta = dsqrt(dabs(D(I_X) - D(I_Z))) * rho * V(:, I_X)
                Vtau = dsqrt(dabs(D(I_Y) - D(I_Z))) * rho * V(:, I_Y)

                Vetas(:, n) = Veta
                Vtaus(:, n) = Vtau
            endif

            do j = 1, 4
                ! support point in the reference frame of the ellipsoid
                ! theta = j * pi / 2 - pi / 4 (j = 1, 2, 3, 4)
                ! x_j = eta * rho * cos(theta_j)
                ! y_j = tau * rho * sin(theta_j)

                ! suppport point in the global reference frame
                points(:, j) = position           &
                             + Veta * costheta(j) &
                             + Vtau * sintheta(j)
            enddo
        end function get_ellipsoid_points

        function get_angles(B, volume) result(angles)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: evec(n_dim, n_dim)
            double precision             :: angles(2) ! (/azimuth, polar/)

            evec = get_eigenvectors(B, volume)

            ! azimuthal angle
            angles(I_X) = datan2(evec(I_Y, I_X), evec(I_X, I_X))

            ! polar angle
            angles(I_Y) = dasin(evec(I_Z, I_Z))

        end function get_angles

end module parcel_ellipsoid
