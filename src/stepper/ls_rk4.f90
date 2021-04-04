module ls_rk4
    use options, only : parcel_info
    implicit none

    contains

        ! allocate memory of temporaries
        subroutine ls_rk4_alloc(num)
            integer, intent(in) :: num

            ! TODO

        end subroutine ls_rk4_alloc

        ! deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            ! TODO

        end subroutine ls_rk4_dealloc


        subroutine ls_rk4_step(dt)
            double precision, intent(in) :: dt

            if (parcel_info%is_elliptic) then
                call ls_rk4_elliptic(dt)
            else
                call ls_rk4_non_elliptic(dt)
            endif

        end subroutine ls_rk4_step


        subroutine ls_rk4_elliptic(dt)
            double precision, intent(in) :: dt

            ! TODO

        end subroutine ls_rk4_elliptic


        subroutine ls_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            ! TODO

        end subroutine ls_rk4_non_elliptic

end module ls_rk4
