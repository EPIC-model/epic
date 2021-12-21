! =============================================================================
!     This module is the 3D version and contains all ellipsoid operations.
!     Reference:
!       Dritschel, D., Reinaud, J., & McKiver, W. (2004).
!       The quasi-geostrophic ellipsoidal vortex model.
!       Journal of Fluid Mechanics, 505, 201-223.
!       doi:10.1017/S0022112004008377
! =============================================================================
module parcel_ellipsoid
    use constants, only : fpi   &
                        , fpi4  &
                        , f34   &
                        , zero  &
                        , two   &
                        , three &
                        , five  &
                        , seven &
                        , max_num_parcels
    use jacobi
    implicit none

    double precision, parameter :: rho = dsqrt(two / five)
    double precision, parameter :: f3pi4 = three * fpi4
    double precision, parameter :: f5pi4 = five * fpi4
    double precision, parameter :: f7pi4 = seven * fpi4
    double precision, parameter :: costheta(4) = dcos((/fpi4, f3pi4, f5pi4, f7pi4/))
    double precision, parameter :: sintheta(4) = dsin((/fpi4, f3pi4, f5pi4, f7pi4/))

    double precision :: etas(max_num_parcels), &
                        taus(max_num_parcels)

    private :: rho, f3pi4, f5pi4, f7pi4, costheta, sintheta, get_upper_triangular

    contains

        ! Obtain the parcel shape matrix.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns the upper trinagular matrix
        function get_upper_triangular(B, volume) result(U)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(3, 3)

            U(1, 1) = B(1)
            U(1, 2) = B(2)
            U(1, 3) = B(3)
            U(2, 2) = B(4)
            U(2, 3) = B(5)
            U(3, 3) = get_B33(B, volume)
        end function get_upper_triangular

        ! Obtain all eigenvalues sorted in descending order
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns all eigenvalues (sorted in descending order)
        function get_eigenvalues(B, volume) result(D)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(3, 3)
            double precision             :: D(3)

            U = get_upper_triangular(B, volume)

            call jacobi_eigenvalues(U, D)

        end function get_eigenvalues

        ! Obtain the normalized eigenvectors.
        ! The eigenvector V(:, j) belongs to the j-th
        ! eigenvalue.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns the eigenvectors
        function get_eigenvectors(B, volume) result(V)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: U(3, 3), D(3), V(3, 3)

            U = get_upper_triangular(B, volume)

            call jacobi_diagonalise(U, D, V)

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
            double precision, intent(out) :: D(3), V(3, 3)
            double precision              :: U(3, 3)

            U = get_upper_triangular(B, volume)

            call jacobi_diagonalise(U, D, V)

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

            if (dabs(B(1) * B(4) - B(2) ** 2) <= epsilon(abc)) then
                print *, "Error in get_B33: Division by small number!"
                stop
            endif

            B33 = (abc ** 2 - B(3) * (B(2) * B(5) - B(3) * B(4)) &
                            + B(1) * B(5) ** 2                   &
                            - B(2) * B(3) * B(5))                &
                / (B(1) * B(4) - B(2) ** 2)

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
            double precision, intent(in) :: D(3)
            double precision             :: lam

            lam = dsqrt(D(1) / D(3))
        end function get_aspect_ratio

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] position vector of the parcel
        ! @param[in] volume of the parcel
        ! @param[in] B matrix elements of the parcel
        ! @returns the parcel support points
        function get_ellipsoid_points(position, volume, B, n, l_reuse) result(points)
            double precision,  intent(in) :: position(3)
            double precision,  intent(in) :: volume
            double precision,  intent(in) :: B(5)        ! B11, B12, B13, B22, B23
            integer, optional, intent(in) :: n
            logical, optional, intent(in) :: l_reuse
            double precision              :: eta, tau, D(3), V(3, 3)
            integer                       :: j
            double precision              :: points(3, 4), xy(2)


            if (present(l_reuse)) then
                if (l_reuse) then
                    eta = etas(n)
                    tau = taus(n)
                else
                    call diagonalise(B, volume, D, V)
                    eta = dsqrt(dabs(D(1) - D(3))) * rho
                    tau = dsqrt(dabs(D(2) - D(3))) * rho

                    etas(n) = eta
                    taus(n) = tau
                endif
            else
                ! (/a2, b2, c2/) with a >= b >= c
                ! D = (/a2, b2, c2/)
                call diagonalise(B, volume, D, V)

                eta = dsqrt(dabs(D(1) - D(3))) * rho
                tau = dsqrt(dabs(D(2) - D(3))) * rho

                etas(n) = eta
                taus(n) = tau
            endif

            do j = 1, 4
                ! support point in the reference frame of the ellipsoid
                ! theta = j * pi / 2 - pi / 4 (j = 1, 2, 3, 4)
                ! x_j = eta * rho * cos(theta_j)
                ! y_j = tau * rho * sin(theta_j)
                xy = (/eta * costheta(j), tau * sintheta(j)/)

                ! suppport point in the global reference frame
                points(:, j) = position        &
                             + xy(1) * V(:, 1) &
                             + xy(2) * V(:, 2)
            enddo
        end function get_ellipsoid_points

        function get_angles(B, volume) result(angles)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: evec(3, 3)
            double precision             :: angles(2) ! (/azimuth, polar/)

            evec = get_eigenvectors(B, volume)

            ! azimuthal angle
            angles(1) = datan2(evec(2, 1), evec(1, 1))

            ! polar angle
            angles(2) = dasin(evec(3, 3))

        end function get_angles

end module parcel_ellipsoid
