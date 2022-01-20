! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the Jacobi algorithm according to
!       Rutishauser, H. The Jacobi method for real symmetric matrices.
!       Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223
! =============================================================================
module jacobi
    use constants, only : hundred, one, zero, f12
    implicit none

    private
        integer, parameter :: n = 3 ! we have 3x3 matrices
        double precision, parameter :: atol = 1.0e-15

    public :: jacobi_diagonalise, jacobi_eigenvalues

    contains

        pure subroutine givens(A, D, i, j, c, s, t, tau)
            double precision, intent(in)  :: A(n, n), D(n)
            integer,          intent(in)  :: i, j
            double precision, intent(out) :: c, s, t, tau
            double precision              :: theta, g, h, aij

            aij = A(i, j)
            g = hundred * dabs(aij)
            h = D(j) - D(i)
            if (dabs(h) + g == dabs(h)) then
                t = aij / h
            else
                theta = f12 * h / aij
                t = one / (dabs(theta) + dsqrt(one + theta ** 2))
                if (theta < zero) then
                    t = -t
                endif
            endif

            ! c = cos(theta)
            c = one / dsqrt(one + t ** 2)

            ! s = sin(theta)
            s = t * c

            tau = s / (one + c)
        end subroutine

        !
        ! Jacobi algorithm -- eigenvalues and eigenvectors
        !

        ! Apply Jacobi (or Givens) rotation to make
        ! A(i, j) = 0. A' = P^T * A * P.
        ! @param[inout] A current state of 3x3 matrix
        ! @param[inout] V current state of eigenvectors
        ! @param[in] i the row
        ! @param[in] j first column with i .ne. j
        ! @param[in] k second column with i .ne. k
        pure subroutine apply_rotation_AV(A, D, Z, i, j, V)
            double precision, intent(inout) :: A(n, n), D(n), Z(n)
            integer,          intent(in)    :: i, j
            double precision, intent(inout) :: V(n, n)
            double precision                :: c, s, t, tau
            integer                         :: k
            double precision                :: g, h

            ! compute the rotation angle theta
            ! Reference:    Rutishauser, H. The Jacobi method for real symmetric matrices.
            !               Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223

            call givens(A, D, i, j, c, s, t, tau)

            !
            ! Apply Givens rotation to matrix
            !

            h = t * A(i, j)
            Z(i) = Z(i) - h
            Z(j) = Z(j) + h
            D(i) = D(i) - h
            D(j) = D(j) + h
            A(i, j) = zero

            do k = 1, i-1
                g = A(k, i)
                h = A(k, j)
                A(k, i) = g - s * (h + g * tau)
                A(k, j) = h + s * (g - h * tau)
            enddo

            do k = i+1, j-1
                g = A(i, k)
                h = A(k, j)
                A(i, k) = g - s * (h + g * tau)
                A(k, j) = h + s * (g - h * tau)
            enddo

            do k = j+1, n
                g = A(i, k)
                h = A(j, k)
                A(i, k) = g - s * (h + g * tau)
                A(j, k) = h + s * (g - h * tau)
            enddo

            do k = 1, n
                ! accumulate eigenvector
                g = V(k, i)
                h = V(k, j)
                V(k, i) = c * g - s * h
                V(k, j) = s * g + c * h
            enddo
        end subroutine apply_rotation_AV

        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten and will be the diagonal
        ! matrix storing the eigenvalues, V contains the eigenvectors.
        ! The eigenvector V(:, i) belongs to eigenvalue A(i, i).
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] D eigenvalues in descending order
        ! @param[out] V eigenvector matrix
        pure subroutine jacobi_diagonalise(A, D, V)
            double precision, intent(inout) :: A(n, n)
            double precision, intent(out)   :: D(n)
            double precision, intent(out)   :: V(n, n)
            double precision                :: B(n), Z(n)
            integer                         :: i, j
            double precision                :: sm

            ! initialise eigenvector matrix to identity matrix
            do j = 1, n
                D(j) = A(j, j)
                B(j) = D(j)
                Z(j) = zero
                do i = 1, n
                    V(i, j) = zero
                    if (i == j) then
                        V(i, j) = one
                    endif
                enddo
            enddo


            ! sum of off-diagonal entries
            ! sm should convergence to zero
            sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))

            do while (sm > atol)

                do i = 1, n-1
                    do j = i+1, n
                        call apply_rotation_AV(A, D, Z, i, j, V)
                    enddo
                enddo

                B = B + Z
                D = B
                Z = zero
                ! update sum of off-diagonals
                sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))
            enddo

            call sort_descending(D, V)

        end subroutine jacobi_diagonalise

        ! Sort the eigenvalues in descending order.
        ! It sorts the eigenvector matrix accordingly.
        ! @param[inout] D eigenvalues
        ! @param[inout] V eigenvector matrix
        pure subroutine sort_descending(D, V)
            double precision, intent(inout) :: D(n)
            double precision, intent(inout) :: V(n, n)
            double precision                :: teval, tevec(n)

            if (D(2) > D(1)) then
                teval = D(1)
                D(1) = D(2)
                D(2) = teval
                tevec = V(:, 1)
                V(:, 1) = V(:, 2)
                V(:, 2) = tevec
            endif

            if (D(3) > D(2)) then
                teval = D(2)
                D(2) = D(3)
                D(3) = teval
                tevec = V(:, 2)
                V(:, 2) = V(:, 3)
                V(:, 3) = tevec
            endif

            if (D(2) > D(1)) then
                teval = D(1)
                D(1) = D(2)
                D(2) = teval
                tevec = V(:, 1)
                V(:, 1) = V(:, 2)
                V(:, 2) = tevec
            endif
        end subroutine sort_descending

        !
        ! Jacobi algorithm -- eigenvalues only
        !

        pure subroutine apply_rotation(A, D, Z, i, j)
            double precision, intent(inout) :: A(n, n), D(n), Z(n)
            integer,          intent(in)    :: i, j
            double precision                :: c, s, t, tau
            integer                         :: k
            double precision                :: g, h

            ! compute the rotation angle theta
            ! Reference:    Rutishauser, H. The Jacobi method for real symmetric matrices.
            !               Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223

            call givens(A, D, i, j, c, s, t, tau)

            !
            ! Apply Givens rotation to matrix
            !

            h = t * A(i, j)
            Z(i) = Z(i) - h
            Z(j) = Z(j) + h
            D(i) = D(i) - h
            D(j) = D(j) + h
            A(i, j) = zero

            do k = 1, i-1
                g = A(k, i)
                h = A(k, j)
                A(k, i) = g - s * (h + g * tau)
                A(k, j) = h + s * (g - h * tau)
            enddo

            do k = i+1, j-1
                g = A(i, k)
                h = A(k, j)
                A(i, k) = g - s * (h + g * tau)
                A(k, j) = h + s * (g - h * tau)
            enddo

            do k = j+1, n
                g = A(i, k)
                h = A(j, k)
                A(i, k) = g - s * (h + g * tau)
                A(j, k) = h + s * (g - h * tau)
            enddo
        end subroutine apply_rotation


        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten, the diagonal
        ! entries will be storing the eigenvalues
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] D eigenvalues in descending order
        pure subroutine jacobi_eigenvalues(A, D)
            double precision, intent(inout) :: A(n, n)
            double precision, intent(out)   :: D(n)
            double precision                :: B(n), Z(n)
            integer                         :: i, j
            double precision                :: sm

            ! initialise
            do j = 1, n
                D(j) = A(j, j)
                B(j) = D(j)
                Z(j) = zero
            enddo

            ! sum of off-diagonal entries
            ! sm should convergence to zero
            sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))

            do while (sm > atol)

                do i = 1, n-1
                    do j = i+1, n
                        call apply_rotation(A, D, Z, i, j)
                    enddo
                enddo

                B = B + Z
                D = B
                Z = zero
                ! update sum of off-diagonals
                sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))
            enddo

            call sort_eigenvalues(D)

        end subroutine jacobi_eigenvalues

        ! Sort the eigenvalues in descending order.
        ! @param[inout] D eigenvalues
        pure subroutine sort_eigenvalues(D)
            double precision, intent(inout) :: D(n)
            double precision                :: teval

            if (D(2) > D(1)) then
                teval = D(1)
                D(1) = D(2)
                D(2) = teval
            endif

            if (D(3) > D(2)) then
                teval = D(2)
                D(2) = D(3)
                D(3) = teval
            endif

            if (D(2) > D(1)) then
                teval = D(1)
                D(1) = D(2)
                D(2) = teval
            endif
        end subroutine sort_eigenvalues

end module jacobi
