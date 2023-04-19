! =============================================================================
!               Timer module implemented according to PMPIC
!                     https://github.com/EPCCed/pmpic
! =============================================================================
module mpi_timer
    use mpi_communicator
    implicit none

    type timer_type
        character(len=32) :: name
        integer           :: handle
        double precision  :: wall_time
        logical           :: running
        integer           :: n_calls
        double precision  :: start_time, end_time
        double precision  :: mean_time
        double precision  :: min_time
        double precision  :: max_time
    end type timer_type

    type(timer_type), allocatable :: timings(:)

    integer :: n_timers = 0
    logical :: l_collected = .false.

    character(7) :: time_unit = 'seconds'
    character(4) :: call_unit = 'one'

    private :: n_timers, time_unit, select_units, resize_timer_by, &
               get_statistics, collect_statistics

    contains

        subroutine resize_timer_by(add)
            integer, intent(in)           :: add
            integer                       :: cur_size
            type(timer_type), allocatable :: tmp(:)

            cur_size = size(timings)

            allocate(tmp(cur_size + add))

            tmp(1:cur_size) = timings

            deallocate(timings)

            ! deallocates tmp
            call move_alloc(tmp, timings)

        end subroutine resize_timer_by

        subroutine register_timer(name, handle)
            character(*), intent(in)  :: name
            integer,      intent(out) :: handle

            if (.not. allocated(timings)) then
                allocate(timings(1))
            endif

            n_timers = n_timers + 1

            if (n_timers > size(timings)) then
                call resize_timer_by(1)
            endif

            handle = n_timers

            timings(handle)%name = name
            timings(handle)%handle = handle
            timings(handle)%wall_time = 0.0d0
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

            timings(handle)%start_time = MPI_Wtime()

        end subroutine start_timer

        subroutine stop_timer(handle)
            integer, intent(in) :: handle

            if (.not. timings(handle)%running) then
                return
            endif

            timings(handle)%running = .false.

            timings(handle)%end_time = MPI_Wtime()

            timings(handle)%wall_time = timings(handle)%wall_time &
                                      + timings(handle)%end_time  &
                                      - timings(handle)%start_time
        end subroutine stop_timer

        subroutine print_timer
            integer          :: n
            double precision :: frac

            call collect_statistics

            if (.not. comm%rank == comm%master) then
                return
            endif

            call select_units

            frac = 100.0d0 / timings(1)%mean_time

            write (*, '("|--------------------------------------------&
                       &----------------------------------------------|")')
            write (*, '("|            Function            &
                        &| #Calls  &
                        &|  % of  &
                        &|  Min time  &
                        &| Mean time  &
                        &|  Max time  |")')

            if (call_unit == 'one') then
                write (*, '("|              name              &
                            &|         &
                            &|  time  &
                            &| in '//time_unit//' &
                            &| in '//time_unit//' &
                            &| in '//time_unit//' |")')
            else
                write (*, '("|              name              &
                            &| in '//call_unit//' &
                            &|  time  &
                            &| in '//time_unit//' &
                            &| in '//time_unit//' &
                            &| in '//time_unit//' |")')
            endif

            write (*, '("|--------------------------------------------&
                       &----------------------------------------------|")')

            do n = 1, n_timers
                write(*,"('|',a,&
                         &'|', i9,&
                         &'| ', f6.2,'%',&
                         &'|', f12.3,&
                         &'|', f12.3,&
                         &'|', f12.3, '|')") &
                timings(n)%name,                &
                timings(n)%n_calls,             &
                frac * timings(n)%mean_time,    &
                timings(n)%min_time,            &
                timings(n)%mean_time,           &
                timings(n)%max_time
            enddo
            write (*, '("|--------------------------------------------&
                       &----------------------------------------------|")')
        end subroutine print_timer

        subroutine write_time_to_csv(fname)
            character(*), intent(in)    :: fname
            integer                     :: n
            double precision            :: frac
            character(len=32)           :: s1, s2, s3, s4, s5
            character(len=len(fname)+4) :: csv
            logical                     :: exists = .false.
            character(len=3)            :: status = 'new'

            call collect_statistics

            if (.not. comm%rank == comm%master) then
                return
            endif

            csv = trim(fname) // '.csv'

            ! check whether file exists
            inquire(file=csv, exist=exists)

            if (exists) then
                status = 'old'
            endif

            call select_units

            frac = 100.0d0 / timings(1)%mean_time

            open(unit=1234, file=csv, status=status)

            if (call_unit == 'one') then
                write(1234, '("function name,&
                              &#calls,&
                              &percentage,&
                              &min time in ' //time_unit//',&
                              &mean time in '//time_unit//',&
                              &max time in ' //time_unit//'")')
            else
                write(1234, '("function name,&
                              &#calls in '//call_unit//',&
                              &percentage,&
                              &min time in ' //time_unit//',&
                              &mean time in '//time_unit//',&
                              &max time in ' //time_unit//'")')
            endif

            do n = 1, n_timers
                write(s1, '(i8)') timings(n)%n_calls
                write(s2, '(f6.2)') frac * timings(n)%mean_time
                write(s3, '(f9.3)') timings(n)%min_time
                write(s4, '(f9.3)') timings(n)%mean_time
                write(s5, '(f9.3)') timings(n)%max_time
                write(1234, "(a,',', a,',', a,',', a,',', a,',', a)") &
                trim(timings(n)%name), &
                trim(adjustl(s1)), &
                trim(adjustl(s2)), &
                trim(adjustl(s3)), &
                trim(adjustl(s4)), &
                trim(adjustl(s5))
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
                call_unit = 'kilo'
                timings(1:n_timers)%n_calls = timings(1:n_timers)%n_calls / 1000
            endif

        end subroutine select_units

        function get_statistics(op) result(buffer)
            type(MPI_Op), intent(in) :: op
            double precision         :: buffer(size(timings))
            integer                  :: buf_size

            buffer = timings(:)%wall_time

            buf_size = size(timings)

            if (comm%rank == comm%master) then
                call MPI_Reduce(MPI_IN_PLACE,           &
                                buffer(1:buf_size),     &
                                buf_size,               &
                                MPI_DOUBLE_PRECISION,   &
                                op,                     &
                                comm%master,            &
                                comm%world,             &
                                comm%err)

            else
                call MPI_Reduce(buffer(1:buf_size),     &
                                buffer(1:buf_size),     &
                                buf_size,               &
                                MPI_DOUBLE_PRECISION,   &
                                op,                     &
                                comm%master,            &
                                comm%world,             &
                                comm%err)
            endif

        end function get_statistics

        subroutine collect_statistics

            if (l_collected) then
                return
            endif

            l_collected = .true.

            ! 3 reductions may be as fast as gather all the data to the root process and to the
            ! calculation there.
            timings(:)%max_time = get_statistics(MPI_MAX)
            timings(:)%min_time = get_statistics(MPI_MIN)
            timings(:)%mean_time = get_statistics(MPI_SUM)

            timings(:)%mean_time = timings(:)%mean_time / dble(comm%size)

        end subroutine collect_statistics

end module mpi_timer
