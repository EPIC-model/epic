module field_mpi
    use mpi_layout
    implicit none

    contains

        subroutine field_halo_fill(data)
            double precision :: data(box%hlo(3):box%hhi(3), &
                                     box%hlo(2):box%hhi(2), &
                                     box%hlo(1):box%hhi(1))

            ! send to left process
            call MPI_Send(data(box%lo(3)

            ! receive from right process


            ! send to right process

            ! receive from left process


            ! send to lower process

            ! receive from upper process

            ! send to upper process

            ! receive from lower process

        end subroutine field_halo_fill

        subroutine field_halo_accumulate(data)

        end subroutine field_halo_accumulate

end module field_mpi
