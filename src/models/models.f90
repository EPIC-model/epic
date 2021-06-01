! =============================================================================
!               This module initializes all parcels and fields.
! =============================================================================
module models
    use parcel_init, only : parcel_default
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
                case ('Robert-uniform')
                    call robert_uniform_init
                case ('Robert-gaussian')
                    call robert_gaussian_init
                case default
                    print *, "Invalid simulation type: '", trim(name), "'"
                    stop
            end select
        end subroutine model_init
end module models
