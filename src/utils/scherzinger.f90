! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the algorithm by Scherzinger and Dohrmann, according to
!       W.M. Scherzinger, C.R. Dohrmann,
!       A robust algorithm for finding the eigenvalues and eigenvectors of
!       3x3 symmetric matrices, Computer Methods in Applied Mechanics and
!       Engineering, Volume 197, Issues 45-48, 2008, Pages 4007-4015,
!       ISSN 0045-7825, https://doi.org/10.1016/j.cma.2008.03.031.
! =============================================================================
module scherzinger
    use constants, only : one, zero, f12, f13, f32, fpi6, three, pi, two, four
    use linalg, only : cross, determinant, trace, signum
    implicit none

    private
        double precision, parameter :: atol = 1.0e-12
        double precision            :: r(3, 3)

    public :: scherzinger_diagonalise, scherzinger_eigenvalues

    contains

        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten, V contains the eigenvectors.
        ! The eigenvector V(:, i) belongs to eigenvalue D(i).
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] D eigenvalues in descending order
        ! @param[out] V eigenvector matrix
        subroutine scherzinger_diagonalise(A, D, V)
            double precision, intent(inout) :: A(3, 3)
            double precision, intent(out)   :: D(3)
            double precision, intent(out)   :: V(3, 3)
            double precision                :: tmp
            double precision                :: rlen(3)
            integer                         :: i

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

            call scherzinger_eigenvalues(A, D)

            V(:, 1) = cross(r(:, 1), r(:, 2))

            ! u1 and u2
            rlen = matmul(A, r(:, 1))
            r(:, 1) = rlen
            rlen = matmul(A, r(:, 2))
            r(:, 2) = rlen

            rlen = norm2(r, dim=1)

            i = maxloc(rlen(1:2), dim=1)

            V(:, 3) = r(:, i) / rlen(i)


            V(:, 2) = cross(V(:, 3), V(:, 1))
            V(:, 3) = cross(V(:, 1), V(:, 2))


            call sort_descending(D, V)

        end subroutine scherzinger_diagonalise

        ! Sort the eigenvalues in descending order.
        ! It sorts the eigenvector matrix accordingly.
        ! @param[inout] D eigenvalues
        ! @param[inout] V eigenvector matrix
        pure subroutine sort_descending(D, V)
            double precision, intent(inout) :: D(3)
            double precision, intent(inout) :: V(3, 3)
            double precision                :: teval, tevec(3)

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
        subroutine scherzinger_eigenvalues(A, D)
            double precision, intent(inout) :: A(3, 3)
            double precision, intent(out)   :: D(3)
            double precision                :: AA(3, 3)
            double precision                :: tr, j2, j3, alpha, eta(3)
            double precision                :: tmp, smp
            double precision                :: s1As1, s1As2, s2As1, s2As2, rlen(3)
            integer                         :: i

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
            alpha = f13 * dacos(f12 * j3 * (three / j2) ** f32)

            if (alpha > fpi6) then
                alpha = alpha + two * f13 * pi
            endif

            eta(1) = two * dsqrt(f13 * j2) * dcos(alpha)

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
            else
                ! i == 3
                rlen = r(:, 3) / rlen(3)

                r(:, 3) = r(:, 1)
                r(:, 1) = rlen
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

!             V(:, 1) = cross(r(:, 1), r(:, 2))

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



            D = eta + f13 * tr

            call sort_eigenvalues(D)

        end subroutine scherzinger_eigenvalues

        ! Sort the eigenvalues in descending order.
        ! @param[inout] D eigenvalues
        pure subroutine sort_eigenvalues(D)
            double precision, intent(inout) :: D(3)
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
