! =============================================================================
!                       Test Jacobi diagonalisation
!
!               This unit test checks the eigenvalue solver for
!               spherical parcels. In the spherical case, the eigenvectors
!               can point in any direction (but they are orthogonal to each
!               other)
! =============================================================================
program test_jacobi_3
    use unit_test
    use constants, only : zero, ten
    use jacobi
    implicit none

    double precision     :: error
    integer              :: n, k, sk
    double precision     :: V(3, 3)
    double precision     :: B(3, 3)
    double precision     :: A(3, 3), eval ! reference
    integer, allocatable :: seed(:)

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)


    error = zero

    do n = 1, 10000

        call create_symmetric_matrix

        call jacobi_diagonalise(B, V)

        ! check eigenvalues
        error = max(error, dabs(eval - B(1, 1)) &
                         + dabs(eval - B(2, 2)) &
                         + dabs(eval - B(3, 3)))

        ! check eigenvectors
        B = matmul(matmul(transpose(V), B), V)

        B = dabs(A - B)

        error = max(error, sum(B))

        ! check if orthogonal
        error = max(error, dabs(dot_product(V(:, 1), V(:, 2))))
        error = max(error, dabs(dot_product(V(:, 1), V(:, 3))))
        error = max(error, dabs(dot_product(V(:, 2), V(:, 3))))
    enddo

    call print_result_dp('Test Jacobi 3', error, atol=dble(1.0e-11))

    contains

        subroutine create_symmetric_matrix
            integer          :: i, j
            double precision :: rndm(3, 3), val
            logical          :: valid
            double precision :: Q(3, 3)

            valid = .false.

            do while (.not. valid)
                call random_number(eval)
                eval = 10 * eval

                ! Create orthogonal vectors using Gram-Schmidt orthogonalization
                ! (https://de.wikipedia.org/wiki/Gram-Schmidtsches_Orthogonalisierungsverfahren)
                do i = 1, 3
                    do j = 1, 3
                        call random_number(val)
                        rndm(i, j) = 10.0d0 * val
                    enddo
                enddo

                Q = 0.0d0

                do j = 1, 3
                    Q(:, j) = rndm(:, j)
                    do i = 1, j - 1
                        Q(:, j) = Q(:, j) &
                                - Q(:, i) * dot_product(Q(:, i), rndm(:, j)) &
                                        / dot_product(Q(:, i), Q(:, i))
                    enddo
                    Q(:, j) = Q(:, j) / norm2(Q(:, j))
                enddo

                ! Check orthogonalization error:
                val = dabs(dot_product(Q(:, 1), Q(:, 2))) &
                    + dabs(dot_product(Q(:, 1), Q(:, 3))) &
                    + dabs(dot_product(Q(:, 2), Q(:, 3)))

                valid = (val < 1.0e-13)
            enddo

            ! Create symmetric B matrix
            B = zero
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        B(i, j) = B(i, j) + Q(i, k) * eval * Q(j, k)
                    enddo
                enddo
            enddo

            A = B
        end subroutine create_symmetric_matrix

end program test_jacobi_3
