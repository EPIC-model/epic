! =============================================================================
!                    Test MPI field halo fill for vector fields
!
!   This unit test checks filling the field halo. Each grid point
!   (including halo) is assigned the value of the rank number + 1. After
!   the call "field_halo_fill", the halo values should be assigned
!   the value of the neighbour rank + 1 owning the halo grid point as
!   interior grid point.
! =============================================================================
program test_field_halo_fill_vector
    use constants, only : zero, one
    use unit_test
    use mpi_environment
    use mpi_layout
    use field_mpi
    implicit none

    integer, parameter            :: nx = 8, ny = 10, nz = 4
    double precision, parameter   :: lower(3) = (/zero, zero, zero/)
    double precision, parameter   :: extent(3) = (/one, one, one/)
    double precision, allocatable :: values(:, :, :, :)
    logical                       :: passed = .true.
    double precision              :: diff, ref, factor
    integer                       :: nc, ncomp

    call mpi_env_initialise

    passed = (world%err == 0)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    diff = zero

    do ncomp = 1, 5

        allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), ncomp))

        do nc = 1, ncomp
            values(:, :, :, nc) = dble(10 ** (nc-1)) * dble(cart%rank + 1)
        enddo

        call field_halo_fill_vector(values, l_alloc=.true.)

        !
        ! check results
        !
        do nc = 1, ncomp
            factor = dble(10 ** (nc-1))

            ! check west halo
            ref = factor * dble(neighbours(MPI_WEST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%hlo(1), nc)))

            ! check the east halo
            ref = factor * dble(neighbours(MPI_EAST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),             &
                                               box%lo(2):box%hi(2),             &
                                               box%hhi(1)-1:box%hhi(1), nc)))

            ! check the south halo
            ref = factor * dble(neighbours(MPI_SOUTH)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),        &
                                               box%hlo(2),                 &
                                               box%lo(1):box%hi(1), nc)))

            ! check the north halo
            ref = factor * dble(neighbours(MPI_NORTH)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),        &
                                               box%hhi(2)-1:box%hhi(2),    &
                                               box%lo(1):box%hi(1), nc)))

            ! check the southwest halo
            ref = factor * dble(neighbours(MPI_SOUTHWEST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),  &
                                               box%hlo(2),           &
                                               box%hlo(1), nc)))

            ! check the northwest halo
            ref = factor * dble(neighbours(MPI_NORTHWEST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),       &
                                               box%hhi(2)-1:box%hhi(2),   &
                                               box%hlo(1), nc)))

            ! check the northeast halo
            ref = factor * dble(neighbours(MPI_NORTHEAST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),            &
                                               box%hhi(2)-1:box%hhi(2),        &
                                               box%hhi(1)-1:box%hhi(1), nc)))

            ! check the southeast halo
            ref = factor * dble(neighbours(MPI_SOUTHEAST)%rank + 1)
            diff = diff + sum(abs(ref - values(box%lo(3):box%hi(3),            &
                                               box%hlo(2),                     &
                                               box%hhi(1)-1:box%hhi(1), nc)))

        enddo

        deallocate(values)
    enddo

    passed = (passed .and. (diff == zero))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI vector field halo fill', passed)
    endif

end program test_field_halo_fill_vector
