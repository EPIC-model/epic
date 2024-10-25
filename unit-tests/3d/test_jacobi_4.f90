! =============================================================================
!                       Test Jacobi diagonalisation
!
!               This unit test checks the eigenvalue solver for
!               for a matrix of the form
!                   /0  0  1\
!                   |0  0  0|
!                   \1  0  0/
!               The eigenvalues of above matrix are +1, 0 and -1 the
!               corresponding (normalised) eigenvectors:
!                   v1 = ( 1/sqrt(2),  0,  1/sqrt(2))
!                   v2 = ( 0,          1,  0)
!                   v3 = (-1/sqrt(2),  0,  1/sqrt(2))
! =============================================================================
program test_jacobi_4
    use unit_test
    use constants, only : zero, one, two
    use jacobi
    implicit none

    double precision     :: error, V(3, 3), B(3, 3), D(3)
    double precision, parameter :: fsqrt2 = one / sqrt(two)
    error = zero

    B = zero
    B(1, 3) = one
    B(3, 1) = one

    call jacobi_diagonalise(B, D, V)

    ! check eigenvalues
    error = max(error, abs( one - D(1)) &
                     + abs(zero - D(2)) &
                     + abs(-one - D(3)))


    error = max(error,                                                  &
                maxval(abs((/fsqrt2, zero, fsqrt2 /) - V(:, 1))),      &
                maxval(abs((/zero,    one,   zero /) - V(:, 2))),      &
                maxval(abs((/fsqrt2, zero, -fsqrt2/) - V(:, 3))))

    call print_result_dp('Test Jacobi 4', error, atol=dble(1.0e-15))

end program test_jacobi_4
