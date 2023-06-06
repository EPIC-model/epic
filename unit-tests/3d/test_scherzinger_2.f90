! =============================================================================
!               Test Scherzinger and Dohrmann diagonalisation
!
!       This unit test checks compares the Jacobi diagonalisation with
!       NumPy (numpy.linalg.eigh) (version 1.20.2)
! =============================================================================
program test_scherzinger_2
    use unit_test
    use constants, only : zero
    use scherzinger
    implicit none

    double precision     :: error
    integer              :: n, k, num
    double precision     :: V(3, 3)
    double precision     :: B(3, 3), D(3)
    double precision     :: Q(3, 3), EV(3) ! reference

    error = zero

    call fopen

    do n = 1, num

        call read_next_line

        call scherzinger_diagonalise(B, D, V)

        ! check eigenvalues
        error = max(error, dabs(EV(1) - D(1)) &
                         + dabs(EV(2) - D(2)) &
                         + dabs(EV(3) - D(3)))


        ! check eigenvectors
        do k = 1, 3
            if (Q(1, k) / V(1, k) < zero) then
                ! opposite sign
                Q(:, k) = - Q(:, k)
            endif
            error = max(error, sum(dabs(Q(:, k) - V(:, k))))
        enddo
    enddo

    call fclose

    call print_result_dp('Test Scherzinger and Dohrmann 2', error, atol=dble(1.0e-12))

    contains

        subroutine fopen
            call open_file('../share/B.asc', 1234)
            read(1234, *) num

            call open_file('../share/V.asc', 1235)
            read(1235, *) n

            if (num .ne. n) then
                ! we do not perform this test
                stop
            endif

            call open_file('../share/D.asc', 1236)
            read(1236, *) n

            if (num .ne. n) then
                ! we do not perform this test
                stop
            endif
        end subroutine fopen

        subroutine open_file(fname, unit)
            character(*), intent(in) :: fname
            integer,      intent(in) :: unit
            logical                  :: exists = .false.

            inquire(file=fname, exist=exists)
            if (.not. exists) then
                call fclose
                ! we do not perform this test
                stop
            endif
            open(unit=unit, file=fname)
        end subroutine open_file

        subroutine fclose
            integer :: unit
            ! 20 October 2021
            ! https://stackoverflow.com/questions/31941681/check-whether-file-has-been-opened-already

            inquire(file='../share/B.asc', number=unit)
            if (unit == 1234) then
                close(1234)
            endif

            inquire(file='../share/V.asc', number=unit)
            if (unit == 1235) then
                close(1235)
            endif

            inquire(file='../share/D.asc', number=unit)
            if (unit == 1236) then
                close(1236)
            endif
        end subroutine fclose

        subroutine read_next_line
            ! B matrix
            read(1234, *) B(1, 1), B(1, 2), B(1, 3), B(2, 2), B(2, 3), B(3, 3)
            B(2, 1) = B(1, 2)
            B(3, 1) = B(1, 3)
            B(3, 2) = B(2, 3)

            ! eigenvectors
            read(1235, *) Q(:, 1), Q(:, 2), Q(:, 3)

            ! eigenvalues
            read(1236, *) EV(1), EV(2), EV(3)
        end subroutine read_next_line

end program test_scherzinger_2
