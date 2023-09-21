module ape_density
    implicit none

    ! If implemented set to .true.
    logical, protected :: l_ape_density = .false.

    public :: ape_den, l_ape_density

    contains

        elemental function ape_den(b, z) result(a)
            double precision, intent(in) :: b       ! buoyancy value
            double precision, intent(in) :: z       ! height
            double precision             :: a       ! APE density

            ! dummy line to avoid compiler warning of 'unused variables'
            a = b + z

        end function ape_den

end module ape_density
