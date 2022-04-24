module mpi_tags
    integer, parameter :: HALO_WEST_TAG      = 1000,  &
                          HALO_EAST_TAG      = 1001,  &
                          HALO_SOUTH_TAG     = 1002,  &
                          HALO_NORTH_TAG     = 1003,  &
                          HALO_SOUTHWEST_TAG = 1004,  &
                          HALO_NORTHWEST_TAG = 1005,  &
                          HALO_NORTHEAST_TAG = 1006,  &
                          HALO_SOUTHEAST_TAG = 1007

    integer, parameter :: SEND_WEST_TAG      = 2000,  &
                          SEND_EAST_TAG      = 2001,  &
                          SEND_SOUTH_TAG     = 2002,  &
                          SEND_NORTH_TAG     = 2003,  &
                          SEND_SOUTHWEST_TAG = 2004,  &
                          SEND_NORTHWEST_TAG = 2005,  &
                          SEND_NORTHEAST_TAG = 2006,  &
                          SEND_SOUTHEAST_TAG = 2007

    integer, parameter :: RECV_WEST_TAG      = 3000,  &
                          RECV_EAST_TAG      = 3001,  &
                          RECV_SOUTH_TAG     = 3002,  &
                          RECV_NORTH_TAG     = 3003,  &
                          RECV_SOUTHWEST_TAG = 3004,  &
                          RECV_NORTHWEST_TAG = 3005,  &
                          RECV_NORTHEAST_TAG = 3006,  &
                          RECV_SOUTHEAST_TAG = 3007

end module mpi_tags
