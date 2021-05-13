! =============================================================================
!               This module initializes all parcels and fields.
! =============================================================================
module model
    use parcel_init, only : parcel_default
    use fields, only : field_default
    use taylorgreen
    implicit none

    contains

        subroutine model_init(name)
            character(*), intent(in) :: name

            call parcel_default

            call field_default

            select case (name)
                case ('TaylorGreen')
                    call taylorgreen_init
                case default
                    print *, "Invalid simulation type: '", name, "'"
                    stop
            end select
        end subroutine model_init
end module model
