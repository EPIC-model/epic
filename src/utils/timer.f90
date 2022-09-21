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

    type(timer_type) :: timings(21)

    integer :: n_timers = 0

    character(7) :: time_unit = 'seconds'
    character(4) :: call_unit = 'one'

    private :: n_timers, time_unit, select_units

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

            call select_units

            frac = hundred / timings(1)%wall_time

            write (*, '("|---------------------------------------------------------------|")')
            write (*, '("|            Function            | #Calls | Total time |  % of  |")')

            if (call_unit == 'one') then
                write (*, '("|              name              |        | in '//time_unit//' |  time  |")')
            else
                write (*, '("|              name              | in '//call_unit// &
                            '| in '//time_unit//' |  time  |")')
            endif

            write (*, '("|---------------------------------------------------------------|")')

            do n = 1, n_timers
                write(*,"('|',a,'|', i8,'|', f12.3, '| ', f6.2,'%','|')") &
                timings(n)%name, timings(n)%n_calls, timings(n)%wall_time, frac * timings(n)%wall_time
            enddo
            write (*, '("|---------------------------------------------------------------|")')
        end subroutine print_timer

        subroutine write_time_to_csv(fname)
            character(*), intent(in)    :: fname
            integer                     :: n
            double precision            :: frac
            character(len=32)           :: s1, s2, s3
            character(len=len(fname)+4) :: csv
            logical                     :: exists = .false.
            character(len=3)            :: status = 'new'

            csv = trim(fname) // '.csv'

            ! check whether file exists
            inquire(file=csv, exist=exists)

            if (exists) then
                status = 'old'
            endif

            call select_units

            frac = hundred / timings(1)%wall_time

            open(unit=1234, file=csv, status=status)

            if (call_unit == 'one') then
                write(1234, '("function name,#calls,total time in ' //time_unit//',percentage")')
            else
                write(1234, '("function name,#calls in '//call_unit//',total time in '//time_unit//',percentage")')
            endif

            do n = 1, n_timers
                write(s1, '(i8)') timings(n)%n_calls
                write(s2, '(f9.3)') timings(n)%wall_time
                write(s3, '(f6.2)') frac * timings(n)%wall_time
                write(1234, "(a,',', a,',', a,',', a)") &
                trim(timings(n)%name), trim(adjustl(s1)), trim(adjustl(s2)), trim(adjustl(s3))
            enddo
            close(1234)

        end subroutine write_time_to_csv

        subroutine select_units
            double precision :: max_wall_time
            integer          :: max_n_calls

            ! go from seconds to minutes (24 hours = 86400 seconds)
            max_wall_time = maxval(timings(1:n_timers)%wall_time)
            if ((max_wall_time > 86400.0d0) .and. (time_unit == 'seconds')) then
                time_unit = 'minutes'
                timings(1:n_timers)%wall_time = timings(1:n_timers)%wall_time / 60.0d0
            endif

            ! go from minutes to hours (7 days = 10080 minutes)
            max_wall_time = maxval(timings(:)%wall_time)
            if ((max_wall_time > 10080.0d0) .and. (time_unit == 'minutes')) then
                time_unit = 'hours'
                timings(1:n_timers)%wall_time = timings(1:n_timers)%wall_time / 60.0d0
            endif

            ! go from hours to days (1 year = 8760 hours)
            max_wall_time = maxval(timings(:)%wall_time)
            if ((max_wall_time > 8760.d0) .and. (time_unit == 'hours')) then
                time_unit = 'days'
                timings(1:n_timers)%wall_time = timings(1:n_timers)%wall_time / 24.0d0
            endif

            max_n_calls = maxval(timings(1:n_timers)%n_calls)
            if ((max_n_calls > 1e6) .and. (call_unit == 'one')) then
                call_unit = 'mega'
                timings(1:n_timers)%n_calls = timings(1:n_timers)%n_calls / 1000
            endif

        end subroutine select_units

end module timer
