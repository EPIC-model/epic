! =============================================================================
!               This module initializes all parcels and fields.
! =============================================================================
module models
    use parcel_init, only : parcel_default, parcel_init_from_grids
    use fields, only : field_default
    use taylorgreen
    use straka
    use robert
    implicit none

    contains

        subroutine model_init(name)
            character(*), intent(in) :: name

            call parcel_default

            call field_default

            select case (trim(name))
                case ('TaylorGreen')
                    call taylorgreen_init
                case ('Straka')
                    call straka_init
                case ('Robert')
                    call robert_init
                case ('FromGrids')
                    call parcel_init_from_grids
                case default
                    print *, "Invalid simulation type: '", trim(name), "'"
                    stop
            end select
        end subroutine model_init
end module models
