! =============================================================================
!           Module to compute 1D derivatives of Fourier series
! =============================================================================
module deriv1d
    use constants, only : pi
    implicit none

    contains

        subroutine init_deriv(n, l, k)
            integer,          intent(in)  :: n
            double precision, intent(in)  :: l
            double precision, intent(out) :: k(n)
            double precision              :: sc
            integer                       :: i

            ! define wavenumber array:
            if (l .ne. 0.d0) then
                ! define wavenumbers:
                sc = pi / l
                do i = 1, n
                    k(i) = sc * dble(i)
                enddo
            else
                ! catastrophic end to run if wave number definition fails:
                write(*,*) '**********************************************'
                write(*,*) ' Wavenumber array definition not possible.'
                write(*,*) ' Domain length equal to zero not allowed.'
                write(*,*) ' STOPPING...'
                write(*,*) '**********************************************'
                stop
            endif
        end subroutine


        ! Calculates m derivatives of length n for the array var(m,n)
        ! and stores the result in der(m,n).
        ! *** Both var and der are assumed to be Fourier ***
        !   *** coeffiecients in their second argument ***
        subroutine deriv(m, n, k, var, der)
            integer,          intent(in)  :: m, n
            double precision, intent(in)  :: k(n), var(m,n)
            double precision, intent(out) :: der(m,n)
            integer                       :: nw, np2, i, j, ic
            double precision              :: ktmp

            nw = n / 2
            np2 = n + 2

            ! Carry out differentiation by wavenumber multiplication:
            do j=1,m
                der(j, 1) = 0.d0
            enddo

            do i = 2, n - nw
                ktmp = k(2*(i-1))
                ic = np2 - i
                do j = 1, m
                    der(j, i) = -ktmp * var(j,ic)
                    der(j,ic) =  ktmp * var(j ,i)
                enddo
            enddo

            if (mod(n,2) .eq. 0) then
                ic=nw+1
                do j = 1, m
                    der(j,ic) = 0.d0
                enddo
            endif
        end subroutine


        ! Calculates m derivatives of length n for the array var(m,n)
        ! and stores the result in der(m,n), where var is represented by a
        ! cosine series. Note der is therefore a sine series.
        ! *** Both var and der are assumed to be Fourier ***
        !   *** coeffiecients in their second argument ***
        subroutine deriv_c(m, n, k, var, der)
            integer,          intent(in)  :: m, n
            double precision, intent(in)  :: k(n), var(m,0:n)
            double precision, intent(out) :: der(m,n)
            integer                       :: i, j
            double precision              :: ktmp

            ! Carry out differentiation by wavenumber multiplication:
            do i = 1, n-1
                ktmp = -k(i)
                do j = 1, m
                    der(j,i) = ktmp * var(j,i)
                enddo
            enddo

            do j = 1, m
                der(j, n) = 0.d0
            enddo
        end subroutine


        ! Calculates m derivatives of length n for the array var(m,n)
        ! and stores the result in der(m,n), where var is represented by a
        ! sine series. Note der is therefore a cosine series.
        ! *** Both var and der are assumed to be Fourier ***
        !   *** coeffiecients in their second argument ***
        subroutine deriv_s(m, n, k, var, der)
            integer,          intent(in)  :: m,n
            double precision, intent(in)  :: k(n), var(m,n)
            double precision, intent(out) :: der(m, 0:n)
            integer                       :: i, j
            double precision              :: ktmp

            ! Carry out differentiation by wavenumber multiplication:
            do i = 1, n-1
                ktmp = k(i)
                do j = 1, m
                    der(j,i) = ktmp * var(j,i)
                enddo
            enddo

            do j = 1, m
                der(j, 0) = 0.d0
                der(j, n) = 0.d0
            enddo
        end subroutine
end module
