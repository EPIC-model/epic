module mpi_tags
    integer, parameter :: HALO_WEST_TAG      = 1000,  &
                          HALO_EAST_TAG      = 1001,  &
                          HALO_SOUTH_TAG     = 1002,  &
                          HALO_NORTH_TAG     = 1003,  &
                          HALO_SOUTHWEST_TAG = 1004,  &
                          HALO_NORTHWEST_TAG = 1005,  &
                          HALO_NORTHEAST_TAG = 1006,  &
                          HALO_SOUTHEAST_TAG = 1007
end module mpi_tags
