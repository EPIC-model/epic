program test_merging_parcels
    use constants, only : pi, zero, one, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_merge_serial
    use parcel_netcdf
    use mpi_communicator
    use mpi_layout
    use mpi_timer
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    call mpi_comm_initialise

    if (comm%size > 1) then
        call mpi_exit_on_error("This program only runs in serial.")
    endif

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)
    call register_timer('parcel I/O', parcel_io_timer)

    parcel%lambda_max = five
    parcel%min_vratio = 8.0d0

    call update_parameters

    call mpi_layout_init(lower, extent, nx, ny, nz)

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


    call create_netcdf_parcel_file('serial_final', .true., .false.)
    call write_netcdf_parcels(t = 0.0d0)

    call mpi_comm_finalise

end program test_merging_parcels
