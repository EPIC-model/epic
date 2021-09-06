! =============================================================================
!               Timer module implemented according to PMPIC
!                     https://github.com/EPCCed/pmpic
! =============================================================================
module timer
    use constants, only : zero, hundred
    implicit none

    type timer_type
        character(len=32) :: name
        integer           :: handle
        double precision  :: wall_time
        logical           :: running
        integer           :: n_calls
        integer(8)        :: cnt1, cnt2
        integer(8)        :: cr1, cr2
    end type timer_type

    type(timer_type) :: timings(18)

    integer :: n_timers = 0

    private :: n_timers

    contains
        subroutine register_timer(name, handle)
            character(*), intent(in)  :: name
            integer,      intent(out) :: handle

            n_timers = n_timers + 1

            handle = n_timers

            timings(handle)%name = name
            timings(handle)%handle = handle
            timings(handle)%wall_time = zero
            timings(handle)%running = .false.
            timings(handle)%n_calls = 0
        end subroutine register_timer

        subroutine start_timer(handle)
            integer, intent(in) :: handle

            if (timings(handle)%running) then
                return
            endif

            timings(handle)%n_calls = timings(handle)%n_calls + 1

            timings(handle)%running = .true.

            call system_clock(timings(handle)%cnt1, timings(handle)%cr1)
        end subroutine start_timer

        subroutine stop_timer(handle)
            integer, intent(in) :: handle

            if (.not. timings(handle)%running) then
                return
            endif

            timings(handle)%running = .false.

            call system_clock(timings(handle)%cnt2, timings(handle)%cr2)

            timings(handle)%wall_time = timings(handle)%wall_time &
                                      + dble(timings(handle)%cnt2) / dble(timings(handle)%cr2) &
                                      - dble(timings(handle)%cnt1) / dble(timings(handle)%cr1)
        end subroutine stop_timer

        subroutine print_timer
            integer          :: n
            double precision :: frac

            write (*, '("|-------------------------------------------------------------|")')
            write (*, '("|            Function            | #Calls |   Total  |  % of  |")')
            write (*, '("|              name              |        |   time   |  time  |")')
            write (*, '("|-------------------------------------------------------------|")')

            frac = hundred / timings(1)%wall_time

            do n = 1, n_timers
                write(*,"('|',a,'|', i8,'|', f9.3,'s', '| ', f6.2,'%','|')") &
                timings(n)%name, timings(n)%n_calls, timings(n)%wall_time, frac * timings(n)%wall_time
            enddo
            write (*, '("|-------------------------------------------------------------|")')
        end subroutine print_timer

        subroutine write_time_to_csv(fname)
            character(*), intent(in)    :: fname
            integer                     :: n
            double precision            :: frac
            character(len=32)           :: s1, s2, s3
            character(len=len(fname)+4) :: csv
            logical                     :: exists = .false.
            character(len=3)            :: status = 'new'

            frac = hundred / timings(1)%wall_time

            csv = trim(fname) // '.csv'

            ! check whether file exists
            inquire(file=csv, exist=exists)

            if (exists) then
                status = 'old'
            endif

            open(unit=1234, file=csv, status=status)

            write(1234, '("function name,#calls,total time,percentage")')
            do n = 1, n_timers
                write(s1, '(i8)') timings(n)%n_calls
                write(s2, '(f9.3)') timings(n)%wall_time
                write(s3, '(f6.2)') frac * timings(n)%wall_time
                write(1234, "(a,',', a,',', a,',', a)") &
                trim(timings(n)%name), trim(adjustl(s1)), trim(adjustl(s2)), trim(adjustl(s3))
            enddo
            close(1234)

        end subroutine write_time_to_csv

end module timer
