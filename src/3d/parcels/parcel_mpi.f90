module parcel_mpi
    use fields, only : get_index
    use field_layout
    implicit none

    private :: locate_parcels,  &
               pack_parcels,    &
               unpack_parcels

    contains

        subroutine parcel_exchange
            integer :: pid(n_parcels)

            ! figure out where parcels should go
            call parcel_locate(pid)

            call pack_parcels

            ! send parcels

            ! receive parcels

            call unpack_parcels

        end subroutine parcel_exchange

        subroutine parcel_locate(pid)
            integer, intent(out) :: pid(:)
            integer              :: i, j, k, n

            do n = 1, n_parcels
                call get_index(parcels%position(:, n), i, j, k)

            enddo

        end subroutine parcel_locate

        subroutine pack_parcels

        end subroutine pack_parcels


        subroutine unpack_parcels

        end subroutine unpack_parcels

end module parcel_mpi
