! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the algorithm by Scherzinger and Dohrmann, according to
!       Rutishauser, H. The Jacobi method for real symmetric matrices.
!       Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223
! =============================================================================
module scherzinger
    use constants, only : hundred, one, zero, f12
    implicit none

    private
        integer, parameter :: n = 3 ! we have 3x3 matrices
        double precision, parameter :: atol = 1.0e-12

    public :: scherzinger_diagonalise, scherzinger_eigenvalues

    contains

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
        pure subroutine scherzinger_diagonalise(A, D, V)
            double precision, intent(inout) :: A(3, 3)
            double precision, intent(out)   :: D(3)
            double precision, intent(out)   :: V(3, 3)
            double precision                :: AA(3, 3)
            double precision                :: tmp, r(3, 3)
            double precision                :: tr, j2, j3, alpha, eta(3)

            ! sum of off-diagonal entries
            tmp = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))

            ! check if a diagonal matrix
            if (tmp < atol) then
                V = zero
                V(1, 1) = one
                V(2, 2) = one
                V(3, 3) = one
                D(1) = A(1, 1)
                D(2) = A(2, 2)
                D(3) = A(3, 3)

                call sort_eigenvalues(D)

                return
            endif

            ! store initial diagonal values
            D(1) = A(1, 1)
            D(2) = A(2, 2)
            D(3) = A(3, 3)

            ! deviatoric matrix:
            tr = trace(A)
            tmp = f13 * tr
            A(1, 1) = A(1, 1) - tmp
            A(2, 2) = A(2, 2) - tmp
            A(3, 3) = A(3, 3) - tmp

            ! invariants j2 and j3:
            AA = matmul(A, A)
            j2 = f12 * trace(AA)
            j3 = determinant(A)

            ! angle:
            alpha = f13 * acos(f12 * j3 * (three / j2) ** f32)

            if (alpha > pi / 6.0d0) then
                alpha = alpha + two * f13 * pi
            endif

            eta(1) = two * sqrt(f13 * j2) * cos(alpha)

            r(:, 1) = (/ A(1, 1) - eta(1), A(2, 1),          A(3, 1)          /)
            r(:, 2) = (/ A(1, 2),          A(2, 2) - eta(1), A(3, 2)          /)
            r(:, 3) = (/ A(1, 3),          A(2, 3),          A(3, 3) - eta(1) /)

            ! length of each r vector
            rlen = norm2(r, dim=1)

            i = maxloc(rlen, dim=1)

            ! s1
            if (i == 1) then
                r(:, 1) = r(:, i) / rlen(i)
            else if (i == 2) then
                rlen = r(:, 2) / rlen(2)

                r(:, 2) = r(:, 1)
                r(:, 1) = rlen
            else if (i == 3) then
                rlen = r(:, 3) / rlen(3)

                r(:, 3) = r(:, 1)
                r(:, 1) = rlen

            else
                print *, "Error"
                stop
            endif

            tmp = dot_product(r(:, 1), r(:, 2))
            smp = dot_product(r(:, 1), r(:, 3))
            r(:, 2) = r(:, 2) - tmp * r(:, 1)
            r(:, 3) = r(:, 3) - smp * r(:, 1)

            ! length of each r vector
            rlen = norm2(r, dim=1)

            i = maxloc(rlen(2:3), dim=1)

            ! s2
            r(:, 2) = r(:, i+1) / rlen(i+1)

            V(:, 1) = cross(r(:, 1), r(:, 2))

            rlen = matmul(A, r(:, 1))
            s1As1 = dot_product(r(:, 1), rlen)
            s2As1 = dot_product(r(:, 2), rlen)

            rlen = matmul(A, r(:, 2))
            s1As2 = dot_product(r(:, 1), rlen)
            s2As2 = dot_product(r(:, 2), rlen)

            tmp = s1As1 - s2As2
            smp = s1As1 + s2As2

            eta(2) = f12 * smp                                                    &
                   - f12 * signum(tmp) * sqrt(tmp ** 2 + four * s1As2 * s2As1)

            eta(3) = smp - eta(2)

            A(1, 1) = D(1) - eta(2)
            A(2, 2) = D(2) - eta(2)
            A(3, 3) = D(3) - eta(2)

            ! u1 and u2
            rlen = matmul(A, r(:, 1))
            r(:, 1) = rlen
            rlen = matmul(A, r(:, 2))
            r(:, 2) = rlen

            rlen = norm2(r, dim=1)

            i = maxloc(rlen(1:2), dim=1)

            V(:, 3) = r(:, i) / rlen(i)


            V(:, 2) = cross(V(:, 3), W(:, 1))
            V(:, 3) = cross(V(:, 1), W(:, 2))

            D = eta + f13 * tr


            call sort_descending(D, V)

        end subroutine scherzinger_diagonalise

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

        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten, the diagonal
        ! entries will be storing the eigenvalues
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] D eigenvalues in descending order
        pure subroutine scherzinger_eigenvalues(A, D)
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

        end subroutine scherzinger_eigenvalues

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

end module scherzinger
