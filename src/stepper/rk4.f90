! =============================================================================
!                       Module for Runge-Kutta methods
! =============================================================================
module rk4
    use classic_rk4
    use ls_rk4
    use constants, only : max_num_parcels
    use options, only : parcel, stepper
    implicit none

    contains

        ! allocate memory of temporaries
        subroutine rk4_alloc(num)
            integer, intent(in) :: num

            if (stepper == 'classic-rk4') then
                call classic_rk4_alloc(num)
            else if (stepper == 'ls-rk4') then
                call ls_rk4_alloc(num)
            else
                print *, "Unknown stepper method '", trim(stepper), "'."
                stop
            endif
        end subroutine rk4_alloc

        ! deallocate memory of temporaries
        subroutine rk4_dealloc
            if (stepper == 'classic-rk4') then
                call classic_rk4_dealloc
            else if (stepper == 'ls-rk4') then
                call ls_rk4_dealloc
            else
                print *, "Unknown stepper method '", trim(stepper), "'."
                stop
            endif
        end subroutine rk4_dealloc


        subroutine rk4_step(dt)
            double precision, intent(in) :: dt

            if (stepper == 'classic-rk4') then
                call classic_rk4_step(dt)
            else if (stepper == 'ls-rk4') then
                call ls_rk4_step(dt)
            else
                print *, "Unknown stepper method '", trim(stepper), "'."
                stop
            endif
        end subroutine rk4_step
end module rk4
