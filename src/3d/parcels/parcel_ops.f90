module parcel_ops
#ifndef NDEBUG
    use parameters, only : lower, upper
#endif
    use parameters, only : extent, extenti
    implicit none

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two zonal coordinates
        ! across periodic edges
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
        ! @returns delx = x1 - x2
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_delx_across_periodic(x1, x2) result (delx)
#else
        elemental function get_delx_across_periodic(x1, x2) result (delx)
#endif
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
#ifndef NDEBUG
            if ((x1 < lower(1)) .or. (x2 < lower(1)) .or. &
                (x1 > upper(1)) .or. (x2 > upper(1))) then
                write(*,*) 'point outside domain was fed into get_delx'
                write(*,*) 'x1, x2, lower(1), upper(1)'
                write(*,*) x1, x2, lower(1), upper(1)
            endif
#endif
            ! works across periodic edge
            delx = delx - extent(1) * dble(nint(delx * extenti(1)))
        end function get_delx_across_periodic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two zonal coordinates
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
        ! @returns delx = x1 - x2
        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx
            delx = x1 - x2
        end function get_delx

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two meridional coordinates
        ! across periodic edges
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_dely_across_periodic(y1, y2) result (dely)
#else
        elemental function get_dely_across_periodic(y1, y2) result (dely)
#endif
            double precision, intent(in) :: y1, y2
            double precision             :: dely

            dely = y1 - y2
#ifndef NDEBUG
            if ((y1 < lower(2)) .or. (y2 < lower(2)) .or. &
                (y1 > upper(2)) .or. (y2 > upper(2))) then
                write(*,*) 'point outside domain was fed into get_dely'
                write(*,*) 'y1, y2, lower(2), upper(2)'
                write(*,*) y1, y2, lower(2), upper(2)
            endif
#endif
            ! works across periodic edge
            dely = dely - extent(2) * dble(nint(dely * extenti(2)))
        end function get_dely_across_periodic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two meridional coordinates
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        elemental function get_dely(y1, y2) result (dely)
            double precision, intent(in) :: y1, y2
            double precision             :: dely
            dely = y1 - y2
        end function get_dely

end module parcel_ops
