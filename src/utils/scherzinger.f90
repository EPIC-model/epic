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
    use constants, only : one, zero, f12, f13, f32, fpi6, three, pi, two, four, f23
    use linalg, only : cross, determinant, trace, signum, dsymm
    implicit none

    private
        double precision, parameter :: atol = 1.0e-12

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
            double precision                :: rlen(3), r(3, 3)
            double precision                :: eta(3)
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

            call evaluate_eigenvalues(A, r, eta, D)

            A(1, 1) = A(1, 1) - eta(2)
            A(2, 2) = A(2, 2) - eta(2)
            A(3, 3) = A(3, 3) - eta(2)

            ! first eigenvector, eq 21:
            V(:, 1) = cross(r(:, 1), r(:, 2))

            ! u1 and u2, eq 30:
            r(:, 1:2) = matmul(A, r(:, 1:2))

            rlen = norm2(r, dim=1)
            i = maxloc(rlen(1:2), dim=1)
            V(:, 3) = r(:, i) / rlen(i)

            ! second and third eigenvector, eq 31:
            V(:, 2) = cross(V(:, 3), V(:, 1))
            V(:, 3) = cross(V(:, 1), V(:, 2))

            if (D(2) > D(1)) then
                tmp = D(1)
                D(1) = D(2)
                D(2) = tmp
                rlen = V(:, 1)
                V(:, 1) = V(:, 2)
                V(:, 2) = rlen
            endif

            if (D(3) > D(2)) then
                tmp = D(2)
                D(2) = D(3)
                D(3) = tmp
                rlen = V(:, 2)
                V(:, 2) = V(:, 3)
                V(:, 3) = rlen
            endif

            if (D(2) > D(1)) then
                tmp = D(1)
                D(1) = D(2)
                D(2) = tmp
                rlen = V(:, 1)
                V(:, 1) = V(:, 2)
                V(:, 2) = rlen
            endif

        end subroutine scherzinger_diagonalise

        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten.
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] D eigenvalues in descending order

        subroutine scherzinger_eigenvalues(A, D)
            double precision, intent(inout) :: A(3, 3)
            double precision, intent(out)   :: D(3)
            double precision                :: eta(3), r(3, 3)
            double precision                :: tmp

            ! sum of off-diagonal entries
            tmp = dabs(A(1, 2)) + dabs(A(1, 3)) + dabs(A(2, 3))

            ! check if a diagonal matrix
            if (tmp < atol) then
                D(1) = A(1, 1)
                D(2) = A(2, 2)
                D(3) = A(3, 3)
            else
                call evaluate_eigenvalues(A, r, eta, D)
            endif

            call sort_eigenvalues(D)

        end subroutine scherzinger_eigenvalues

        subroutine evaluate_eigenvalues(A, r, eta, D)
            double precision, intent(inout) :: A(3, 3)
            double precision, intent(out)   :: eta(3), r(3, 3)
            double precision, intent(out)   :: D(3)
            double precision                :: AA(3, 3)
            double precision                :: tr, j2, j3, alpha
            double precision                :: tmp, smp
            double precision                :: s1As1, s1As2, s2As1, s2As2, rlen(3)
            integer                         :: i

            ! deviatoric matrix (eq 2):
            tr = trace(A)

            tmp = f13 * tr
            A(1, 1) = A(1, 1) - tmp
            A(2, 2) = A(2, 2) - tmp
            A(3, 3) = A(3, 3) - tmp

            ! invariants j2 and j3 (eq 5):
            AA = dsymm(A, A)
            j2 = f12 * trace(AA)
            j3 = determinant(A)

            ! angle (solve eq 7 for alpha):
            alpha = f13 * dacos(f12 * j3 * (three / j2) ** f32)

            if (alpha > fpi6) then
                alpha = alpha + f23 * pi
            endif

            ! eq 6:
            eta(1) = two * dsqrt(f13 * j2) * dcos(alpha)

            ! eq 19:
            r(:, 1) = (/ A(1, 1) - eta(1), A(2, 1),          A(3, 1)          /)
            r(:, 2) = (/ A(1, 2),          A(2, 2) - eta(1), A(3, 2)          /)
            r(:, 3) = (/ A(1, 3),          A(2, 3),          A(3, 3) - eta(1) /)

            ! length of each r vector
            rlen = norm2(r, dim=1)

            i = maxloc(rlen, dim=1)

            ! s1
            rlen = r(:, i) / rlen(i)

            if (i == 2) then
                r(:, 2) = r(:, 1)
            else if (i == 3) then
                r(:, 3) = r(:, 1)
            endif

            ! s1
            r(:, 1) = rlen

            ! t2 and t3 of eq 20:
            tmp = dot_product(r(:, 1), r(:, 2))
            smp = dot_product(r(:, 1), r(:, 3))
            r(:, 2) = r(:, 2) - tmp * r(:, 1)
            r(:, 3) = r(:, 3) - smp * r(:, 1)

            ! length of each r vector
            rlen = norm2(r, dim=1)

            ! s2
            i = maxloc(rlen(2:3), dim=1)
            r(:, 2) = r(:, i+1) / rlen(i+1)

            ! first eigenvector, eq 21:
!             V(:, 1) = cross(r(:, 1), r(:, 2))

            ! eq 22:
            AA(:, 1:2) = matmul(A, r(:, 1:2))
            s1As1 = dot_product(r(:, 1), AA(:, 1))
            s2As1 = dot_product(r(:, 2), AA(:, 1))
            s1As2 = s2As1 !dot_product(r(:, 1), rlen)
            s2As2 = dot_product(r(:, 2), AA(:, 2))


            tmp = s1As1 - s2As2
            smp = s1As1 + s2As2

            ! second eigenvalue, eq 27:
            eta(2) = f12 * smp                                                    &
                   - f12 * signum(tmp) * dsqrt(tmp ** 2 + four * s1As2 * s2As1)

            !  third eigenvalue, eq 27:
            eta(3) = smp - eta(2)

            ! final eigenvalues:
            D = eta + f13 * tr

        end subroutine evaluate_eigenvalues

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
