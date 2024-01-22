module parcel_ops
    use parameters, only : extent, extenti, center, lower, upper
    implicit none

    contains

        ! Obtain the difference between two horizontal coordinates
        ! across periodic edges
        ! @param[in] x1 first horizontal position
        ! @param[in] x2 second horizontal position
        ! @returns delx = x1 - x2
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_delx(x1, x2) result (delx)
#else
        elemental function get_delx(x1, x2) result (delx)
#endif
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
#ifndef NDEBUG
            if ((x1 < lower(1)) .or. (x2 < lower(1)) .or. (x1 > upper(1)) .or. (x2 > upper(1))) then
                write(*,*) 'point outside domain was fed into get_delx'
                write(*,*) 'x1, x2, lower(1), upper(1)'
                write(*,*) x1, x2, lower(1), upper(1)
            endif
#endif
            ! works across periodic edge
            delx = delx - extent(1) * dble(nint(delx * extenti(1)))
        end function get_delx

end module parcel_ops
