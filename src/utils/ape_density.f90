module ape_density
    use constants
    implicit none

    ! If implemented set to .true.
    logical, protected :: l_ape_density = .true.

    public :: ape_den, l_ape_density

    contains

        elemental function ape_den(b, z) result(a)
            double precision, intent(in) :: b       ! buoyancy value
            double precision, intent(in) :: z       ! height
            double precision             :: a       ! APE density
            double precision             :: br

            br = max(b, -twopi)
            br = min(br, twopi)
            a = f18 * (br - four * z) ** 2

        end function ape_den

end module ape_density
