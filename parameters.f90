! global parameter settings
module parameters
    implicit none

    logical :: verbose = .false.

    character(len=32) :: filename = ''

    integer :: h5freq = 1

end module parameters
