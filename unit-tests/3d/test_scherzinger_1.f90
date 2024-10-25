! =============================================================================
!               Test Scherzinger and Dohrmann diagonalisation
!
!               This unit test checks the eigenvalue solver.
! =============================================================================
program test_scherzinger_1
    use unit_test
    use constants, only : zero
    use scherzinger
    implicit none

    double precision     :: error
    integer              :: n, k, sk
    logical              :: valid
    double precision     :: V(3, 3)
    double precision     :: B(3, 3), D(3)
    double precision     :: Q(3, 3), EV(3) ! reference
    integer, allocatable :: seed(:)

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)


    error = zero

    n = 1
    do while (n <= 100000)

        call create_symmetric_matrix

        if (valid) then
            n = n + 1

            call scherzinger_diagonalise(B, D, V)

            ! check eigenvalues
            error = max(error, abs(EV(1) - D(1)) &
                             + abs(EV(2) - D(2)) &
                             + abs(EV(3) - D(3)))


            ! check eigenvectors
            do k = 1, 3
                if (Q(1, k) / V(1, k) < zero) then
                    ! opposite sign
                    Q(:, k) = - Q(:, k)
                endif
                error = max(error, sum(abs(Q(:, k) - V(:, k))))
            enddo
        endif
    enddo

    call print_result_dp('Test Scherzinger and Dohrmann 1', error, atol=dble(1.0e-9))

    contains

        subroutine create_symmetric_matrix
            integer          :: i, j
            double precision :: rndm(3, 3), val
            double precision :: tmp

            ! random eigenvalue
            do j = 1, 3
                call random_number(tmp)
                EV(j) = 2.0d0 * tmp
            enddo

            ! sort
            if (EV(1) < EV(2)) then
                tmp = EV(1)
                EV(1) = EV(2)
                EV(2) = tmp
            endif

            if (EV(2) < EV(3)) then
                tmp = EV(2)
                EV(2) = EV(3)
                EV(3) = tmp
            endif

            if (EV(1) < EV(2)) then
                tmp = EV(1)
                EV(1) = EV(2)
                EV(2) = tmp
            endif


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
            val = abs(dot_product(Q(:, 1), Q(:, 2))) &
                + abs(dot_product(Q(:, 1), Q(:, 3))) &
                + abs(dot_product(Q(:, 2), Q(:, 3)))

            valid = (val < 1.0e-13)
            if (.not. valid) then
                return
            endif

            ! Create symmetric B matrix
            B = 0.0d0
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        B(i, j) = B(i, j) + Q(i, k) * EV(k) * Q(j, k)
                    enddo
                enddo
            enddo
        end subroutine create_symmetric_matrix

end program test_scherzinger_1
