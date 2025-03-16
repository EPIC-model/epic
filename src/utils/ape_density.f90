module ape_density
    use constants
    implicit none

    ! If implemented set to .true.
    logical, protected :: l_ape_density = .false.

    public :: ape_den, l_ape_density

    contains

        elemental function ape_den(b, z) result(a)
            double precision, intent(in) :: b       ! buoyancy value
            double precision, intent(in) :: z       ! height
            double precision             :: a       ! APE density
            !double precision             :: br

            ! dummy line to avoid compiler warning of 'unused variables'
            a = b + z

            ! RT test case
            !br = max(b, -one)
            !br = min(br, one)
            !a = br * dasin(br) + dsqrt(one - br ** 2) - z * br - dcos(z)
            
            ! IW test case
            !a = f18 * (b - four * z) ** 2

        end function ape_den

end module ape_density
