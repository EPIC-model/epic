module field_ops
    use mpi_environment
    use mpi_layout
    use constants, only : f12
    use parameters, only : nz, ncell
    use mpi_utils, only : mpi_check_for_error
    implicit none

    contains

        function get_mean(ff) result(mean)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: mean

            ! (divide by ncell since lower and upper edge weights are halved)
            mean = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                            + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                        + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))) / dble(ncell)

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               mean,                    &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, "in MPI_Allreduce of field_ops::get_mean.")

        end function get_mean


        function get_sum(ff) result(res)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: res

            res = sum(ff(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               res,                     &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, "in MPI_Allreduce of field_ops::get_sum.")

        end function get_sum


        function get_rms(ff) result(rms)
            double precision, intent(in) :: ff(box%hlo(3):box%hhi(3), &
                                               box%hlo(2):box%hhi(2), &
                                               box%hlo(1):box%hhi(1))
            double precision :: rms

            rms = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2  &
                           + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2) &
                       + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)) / dble(ncell)

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               rms,                     &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, "in MPI_Allreduce of field_ops::get_rms.")

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


            call MPI_Allreduce(MPI_IN_PLACE,            &
                               abs_max,                 &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, "in MPI_Allreduce of field_ops::get_abs_max.")

        end function get_abs_max

end module field_ops
