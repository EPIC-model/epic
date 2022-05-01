module field_ops
    use mpi_communicator
    use mpi_layout
    use constants, only : f12
    use parameters, only : nz, ncell
    implicit none

    contains

        function get_mean(ff) result(mean)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: mean

            mean = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                            + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                        + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))) / dble(ncell)

            call MPI_Allreduce(MPI_IN_PLACE, mean, 1, MPI_DOUBLE, MPI_SUM, comm_world, mpi_err)

        end function get_mean


        function get_rms(ff) result(rms)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: rms

            rms = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2  &
                           + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2) &
                       + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)) / dble(ncell)

            call MPI_Allreduce(MPI_IN_PLACE, rms, 1, MPI_DOUBLE, MPI_SUM, comm_world, mpi_err)

            rms = dsqrt(rms)

        end function get_rms

        function get_abs_max(ff) result(abs_max)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: abs_max

            abs_max = maxval(dabs(ff(box%lo(3):box%hi(3),   &
                                     box%lo(2):box%hi(2),   &
                                     box%lo(1):box%hi(1))))


            call MPI_Allreduce(MPI_IN_PLACE, abs_max, 1, MPI_DOUBLE, MPI_MAX, comm_world, mpi_err)

        end function get_abs_max

end module field_ops
