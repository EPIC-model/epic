! =============================================================================
!     This module is the 3D version and contains all ellipsoid operations.
! =============================================================================
module parcel_ellipsoid
    use constants, only : fpi   &
                        , fpi4  &
                        , f34   &
                        , zero  &
                        , two   &
                        , three &
                        , five  &
                        , seven
    use jacobi
    implicit none

    double precision, parameter :: rho = dsqrt(two / five)
    double precision, parameter :: f3pi4 = three * fpi4
    double precision, parameter :: f5pi4 = five * fpi4
    double precision, parameter :: f7pi4 = seven * fpi4
    double precision, parameter :: costheta(4) = dcos((/fpi4, f3pi4, f5pi4, f7pi4/))
    double precision, parameter :: sintheta(4) = dsin((/fpi4, f3pi4, f5pi4, f7pi4/))

    private :: rho

    contains

        ! Obtain the parcel shape matrix.
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns the symmetric 3x3 shape matrix
        function get_symmetric_matrix(B, volume) result(D)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: D(3, 3)

            D(1, 1) = B(1)
            D(1, 2) = B(2)
            D(2, 1) = B(2)
            D(1, 3) = B(3)
            D(3, 1) = B(3)
            D(2, 2) = B(4)
            D(2, 3) = B(5)
            D(3, 2) = B(5)
            D(3, 3) = get_B33(B, volume)
        end function get_symmetric_matrix

        ! Obtain all eigenvalues sorted in descending order
        ! @param[in] B = (B11, B12, B13, B22, B23)
        ! @param[in] volume of the parcel
        ! @returns all eigenvalues (sorted in descending order)
        function get_eigenvalues(B, volume) result(evals)
            double precision, intent(in) :: B(5)
            double precision, intent(in) :: volume
            double precision             :: D(3, 3)
            double precision             :: evals(3)

            D = get_symmetric_matrix(B, volume)

            call jacobi_diagonalise(D)

            evals(1) = D(1, 1)
            evals(2) = D(2, 2)
            evals(3) = D(3, 3)

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
            double precision             :: D(3, 3), V(3, 3)

            D = get_symmetric_matrix(B, volume)

            call jacobi_diagonalise(D, V)

        end function get_eigenvectors

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

            B33 = (abc ** 2 + B(2) * B(5) * (B(5) - two * B(3)) + B(3) ** 2 * B(4)) &
                / (B(1) * B(4) - B(2) ** 2)

        end function get_B33

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] volume of the parcel
        ! @returns abc = 3 * volume / (4 * pi)
        elemental function get_abc(volume) result(abc)
            double precision, intent(in) :: volume
            double precision             :: abc

            abc = f34 * volume * fpi
        end function get_abc

        ! Obtain the aspect ratios a/b of the parcel(s).
        ! @param[in] a2 is the largest eigenvalue
        ! @param[in] b2 is the middle eigenvalue
        ! @param[in] volume of the parcel(s)
        ! @returns (a/b, a/c)
        elemental function get_aspect_ratio(a2, b2) result(lam)
            double precision, intent(in) :: a2, b2
            double precision             :: lam

            lam = dsqrt(a2 / b2)
        end function get_aspect_ratio

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] position vector of the parcel
        ! @param[in] volume of the parcel
        ! @param[in] B matrix elements of the parcel
        ! @returns the parcel support points
        function get_ellipsoid_points(position, volume, B) result(points)
            double precision, intent(in) :: position(3)
            double precision, intent(in) :: volume
            double precision, intent(in) :: B(5)        ! B11, B12, B13, B22, B23
            double precision             :: B33
            double precision             :: eta, tau, evals(3)
            integer                      :: j
            double precision             :: points(4, 3)

            ! (/a2, b2, c2/)
            evals = get_eigenvalues(B, volume)

            eta = dsqrt(evals(2) - evals(1)) * rho
            tau = dsqrt(evals(3) - evals(1)) * rho

            do j = 1, 4
                ! theta = j * pi / 2 - pi / 4 (j = 1, 2, 3, 4)
                ! x_j = eta * rho * cos(theta_j)
                ! y_j = tau * rho * sin(theta_j)
                points(4, :) = (/eta * costheta(1), tau * sintheta(1), zero/)
            enddo

        end function get_ellipsoid_points
end module parcel_ellipsoid
