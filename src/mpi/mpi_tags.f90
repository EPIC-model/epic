module mpi_tags
    integer, parameter :: HALO_WEST_TAG      = 1000,  &
                          HALO_EAST_TAG      = 1001,  &
                          HALO_SOUTH_TAG     = 1002,  &
                          HALO_NORTH_TAG     = 1003,  &
                          HALO_SOUTHWEST_TAG = 1004,  &
                          HALO_NORTHWEST_TAG = 1005,  &
                          HALO_NORTHEAST_TAG = 1006,  &
                          HALO_SOUTHEAST_TAG = 1007

    integer, parameter :: NEIGHBOUR_TAG(8) = (/2000,    &
                                               2001,    &
                                               2002,    &
                                               2003,    &
                                               2004,    &
                                               2005,    &
                                               2006,    &
                                               2007/)

    integer, parameter :: REVERSE_LO_TAG = 3000, &
                          REVERSE_HI_TAG = 3001

end module mpi_tags
