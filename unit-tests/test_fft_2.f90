! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_fft_2
    use unit_test
    use constants, only : pi, twopi, zero, one
    use stafft, only : initfft, forfft
    implicit none

    integer, parameter          :: n = 16
    double precision            :: trig(2 * n)
    integer                     :: factors(5)
    double precision            :: fun(n)
    double precision            :: ref(2 * n)            ! reference solution in Fourier space
    double precision            :: a ! 0 < a < 1
    double precision            :: da = 0.1d0
    double precision            :: err = zero
    double precision, parameter :: atol = dble(1.0e-14)

    ! set up FFT
    call initfft(n, factors, trig)

    ! test function 1
    a = da
    do while (a <= 0.9d0)

        ! initialize function
        fun = func_1(a, n)

        ! compute reference solution
        call dft(n, fun, ref)

        ! forward FFT (Hermitian form) (overwrites fun)
        call forfft(1, n, fun, trig, factors)

        err = max(err, get_max_absolute_error(fun, ref, n))

        a = a + da
    enddo

    ! test function 2
    a = da
    do while (a <= 0.9d0)

        ! initialize function
        fun = func_2(a, n)

        ! compute reference solution
        call dft(n, fun, ref)

        ! forward FFT (Hermitian form) (overwrites fun)
        call forfft(1, n, fun, trig, factors)

        err = max(err, get_max_absolute_error(fun, ref, n))

        a = a + da
    enddo

    ! final check

    call print_result_dp('Test FFT transform', err, atol)

    contains

        ! initialize function in real space from -pi to pi
        ! f(x) = 1 / (1 - a * cos(x + 1))
        ! 0 < a < 1
        function func_1(a, n) result(f1)
            double precision, intent(in) :: a
            integer,          intent(in) :: n
            double precision             :: f1(n)
            double precision             :: dx, xx, phase
            integer                      :: j

            phase = one ! phase offset in radians

            dx = twopi / dble(n)
            do j = 0, n-1
                xx = -pi + dx * dble(j)
                f1(j+1) = one / (one - a * dcos(xx + phase))
            enddo
        end function func_1

        ! initialize function in real space from -pi to pi
        ! f(x) = 1 / (1 - a * sin(x + 1))
        ! 0 < a < 1
        function func_2(a, n) result(f2)
            double precision, intent(in) :: a
            integer,          intent(in) :: n
            double precision             :: f2(n)
            double precision             :: dx, xx, phase
            integer                      :: j

            phase = one ! phase offset in radians

            dx = twopi / dble(n)
            do j = 0, n-1
                xx = -pi + dx * dble(j)
                f2(j+1) = one / (one - a * dsin(xx + phase))
            enddo
        end function func_2

        ! discrete Fourier transform
        subroutine dft(m, f, fhat)
            integer,          intent(in)  :: m
            double precision, intent(in)  :: f(m)       ! function in real space
            double precision, intent(out) :: fhat(2*m)  ! function in Fourier space
            double precision              :: fac, arg
            integer                       :: k, l
            double precision              :: sfac       ! normalization factor

            fhat = 0.d0
            fac = twopi / dble(m)

            do k = 0, m - 1
                do l = 0, m - 1
                    arg = fac * dble(k * l)
                    fhat(k+1)   = fhat(k+1)   + f(l+1) * dcos(arg)
                    fhat(k+1+m) = fhat(k+1+m) - f(l+1) * dsin(arg)
                enddo
            enddo

            sfac = one / dsqrt(dble(m))

            fhat = fhat * sfac

        end subroutine dft


        function get_max_absolute_error(fun, ref, n) result(err)
            double precision, intent(in) :: fun(n), ref(2 * n)
            integer,          intent(in) :: n
            double precision             :: err, err_cos, err_sin
            integer                      :: i, j, k, nw

            nw = n / 2

            ! Cosine coefficients
            err_cos = maxval(abs(ref(1:nw+1) - fun(1:nw+1)))

            ! Sine coefficients
            i = 1
            j = nw - 1
            k = n + 1

            err_sin = maxval(abs(ref(k+i:k+j) - fun(k-i:k-j:-1)))

            err = max(err_cos, err_sin)

        end function get_max_absolute_error

end program test_fft_2
