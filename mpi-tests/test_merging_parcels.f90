program test_merging_parcels
    use constants, only : pi, zero, one, two, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_merge
    use parcel_netcdf
    use mpi_communicator
    use mpi_layout
    use test_utils
    implicit none

    call mpi_comm_initialise

    call register_all_timers

    nx = 32
    ny = 32
    nz = 32
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    parcel%lambda_max = 4.0d0
    parcel%min_vratio = 40.0d0
    parcel%size_factor = 4

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call nearest_win_allocate

    call parcel_alloc(max_num_parcels)

    call read_netcdf_parcels('initial_parcels.nc')

    n_total_parcels = 0

    call MPI_Allreduce(n_parcels,       &
                       n_total_parcels, &
                       1,               &
                       MPI_INTEGER,     &
                       MPI_SUM,         &
                       comm%world,      &
                       comm%err)

    call merge_parcels(parcels)

    n_total_parcels = 0
    call MPI_Allreduce(n_parcels,       &
                       n_total_parcels, &
                       1,               &
                       MPI_INTEGER,     &
                       MPI_SUM,         &
                       comm%world,      &
                       comm%err)

    call create_netcdf_parcel_file('parallel_final', .true., .false.)
    call write_netcdf_parcels(t = 0.0d0)

    call nearest_win_deallocate

    call mpi_comm_finalise

end program test_merging_parcels
