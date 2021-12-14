! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the Jacobi algorithm, see for example
! https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
! =============================================================================
module jacobi
    use constants, only : hundred, one, zero, f12
    implicit none

    private
        integer, parameter :: n = 3 ! we have 3x3 matrices
        double precision, parameter :: atol = 1.0e-15
        logical :: do_evec

    public :: jacobi_diagonalise

    contains

        ! Get the index of the maximum absolute
        ! off-diagonal matrix element.
        ! @param A[in] current state of 3x3 matrix
        ! @param i[in] the row
        ! @param j[in] first column with i .ne. j
        ! @param k[in] second column with i .ne. k
        ! @param l[out] the index of the pivot
        pure function get_pivot(A, i, j, k) result (l)
            double precision, intent(in) :: A(n, n)
            integer,          intent(in) :: i, j, k
            integer                      :: l

            l = k
            if (dabs(A(i, j)) > dabs(A(i, k))) then
                l = j
            endif
        end function get_pivot

        ! Apply Jacobi (or Givens) rotation to make
        ! A(i, j) = 0. A' = P^T * A * P.
        ! @param[inout] A current state of 3x3 matrix
        ! @param[inout] V current state of eigenvectors
        ! @param[in] i the row
        ! @param[in] j first column with i .ne. j
        ! @param[in] k second column with i .ne. k
        subroutine apply_rotation(A, V, i, j)
            double precision, intent(inout)           :: A(n, n)
            double precision, intent(inout), optional :: V(n, n)
            integer,          intent(in)              :: i, j
            double precision                          :: theta, c, s, t, tau
            integer                                   :: k
            double precision                          :: g, h, aij

            ! compute the rotation angle theta
            ! Reference:    Rutishauser, H. The Jacobi method for real symmetric matrices.
            !               Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223
            g = hundred * dabs(A(i, j))
            h = A(j, j) - A(i, i)
            if (dabs(h) + g == dabs(h)) then
                t = A(i, j) / h
            else
                theta = f12 * h / A(i, j)
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

            ! accumulate eigenvector
            if (do_evec) then
                do k = 1, n
                    g = V(i, k)
                    h = V(j, k)
                    V(i, k) = c * g - s * h
                    V(j, k) = s * g + c * h
                enddo
            endif

            !
            ! Apply Givens rotation to matrix
            !

            aij = A(i, j)

            ! update off-diagonal entries
            do k = 1, n
                g = A(k, i)
                h = A(k, j)
                if (.not. (k == i)) then
                    A(k, i) = g - s * (h + tau * g)
                    A(i, k) = A(k, i)
                endif

                if (.not. (k == j)) then
                    A(k, j) = h + s * (g - tau * h)
                    A(j, k) = A(k, j)
                endif
            enddo

            ! update diagonal entries
            A(i, i) = A(i, i) - t * aij
            A(j, j) = A(j, j) + t * aij

            ! set A(i, j) = A(j, i) explicitly to zero
            A(i, j) = zero
            A(j, i) = zero

        end subroutine apply_rotation


        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten and will be the diagonal
        ! matrix storing the eigenvalues, V contains the eigenvectors.
        ! The eigenvector V(:, i) belongs to eigenvalue A(i, i).
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] V eigenvector matrix (optional)
        subroutine jacobi_diagonalise(A, V)
            double precision, intent(inout)           :: A(n, n)
            double precision, intent(out), optional   :: V(n, n)
            integer                                   :: i, j
            double precision                          :: sm

            do_evec = present(V)

            ! initialise eigenvector matrix to identity matrix
            if (do_evec) then
                do j = 1, n
                    do i = 1, n
                        V(i, j) = zero
                        if (i == j) then
                            V(i, j) = one
                        endif
                    enddo
                enddo
            endif


            ! sum of off-diagonal entries
            ! sm should convergence to zero
            sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))

            do while (sm > atol)
                ! first row
                j = get_pivot(A, 1, 2, 3)
                call apply_rotation(A, V, 1, j)

                ! second row
                j = get_pivot(A, 2, 1, 3)
                call apply_rotation(A, V, 2, j)

                ! update sum of off-diagonals
                sm = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))
            enddo

            call sort_descending(A, V)


            if (do_evec) then
                ! eigenvector of i-th eigenvalue
                ! is given by V(:, i)
                V = transpose(V)
            endif

        end subroutine jacobi_diagonalise


        ! sort eigenvalues and eigenvectors
        ! Sort the eigenvalues in descending order.
        ! It sorts the eigenvector matrix accordingly.
        ! @param[inout] D diagonal matrix with eigenvalues
        ! @param[inout] V eigenvector matrix
        subroutine sort_descending(D, V)
            double precision, intent(inout)           :: D(n, n)
            double precision, intent(inout), optional :: V(n, n)
            double precision                          :: teval, tevec(n)

            if (D(2, 2) > D(1, 1)) then
                teval = D(1, 1)
                D(1, 1) = D(2, 2)
                D(2, 2) = teval
                if (do_evec) then
                    tevec = V(1, :)
                    V(1, :) = V(2, :)
                    V(2, :) = tevec
                endif
            endif

            if (D(3, 3) > D(2, 2)) then
                teval = D(2, 2)
                D(2, 2) = D(3, 3)
                D(3, 3) = teval
                if (do_evec) then
                    tevec = V(2, :)
                    V(2, :) = V(3, :)
                    V(3, :) = tevec
                endif
            endif

            if (D(2, 2) > D(1, 1)) then
                teval = D(1, 1)
                D(1, 1) = D(2, 2)
                D(2, 2) = teval
                if (do_evec) then
                    tevec = V(1, :)
                    V(1, :) = V(2, :)
                    V(2, :) = tevec
                endif
            endif
        end subroutine sort_descending

end module jacobi
