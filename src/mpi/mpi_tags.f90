module mpi_tags
    implicit none

    integer, parameter  :: MPI_NONE      = 0, &
                           MPI_NORTH     = 1, &
                           MPI_SOUTH     = 2, &
                           MPI_WEST      = 3, &
                           MPI_EAST      = 4, &
                           MPI_NORTHWEST = 5, &
                           MPI_NORTHEAST = 6, &
                           MPI_SOUTHWEST = 7, &
                           MPI_SOUTHEAST = 8

    integer, parameter :: SEND_NEIGHBOUR_TAG(8) = (/MPI_SOUTH,      &   ! north maps to south
                                                    MPI_NORTH,      &   ! south maps to north
                                                    MPI_EAST,       &   ! west maps to east
                                                    MPI_WEST,       &   ! east maps to west
                                                    MPI_SOUTHEAST,  &   ! northwest maps to southeast
                                                    MPI_SOUTHWEST,  &   ! northeast maps to southwest
                                                    MPI_NORTHEAST,  &   ! southwest maps to northeast
                                                    MPI_NORTHWEST/)

    integer, parameter :: RECV_NEIGHBOUR_TAG(8) = (/MPI_NORTH,      &
                                                    MPI_SOUTH,      &
                                                    MPI_WEST,       &
                                                    MPI_EAST,       &
                                                    MPI_NORTHWEST,  &
                                                    MPI_NORTHEAST,  &
                                                    MPI_SOUTHWEST,  &
                                                    MPI_SOUTHEAST/)

    integer, parameter :: REVERSE_LO_TAG = 3000, &
                          REVERSE_HI_TAG = 3001

end module mpi_tags
