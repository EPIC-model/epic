module field_io
#ifdef ENABLE_HDF5
    use field_hdf5
    use h5_reader
#endif
#ifdef ENABLE_NETCDF
    use field_netcdf
    use netcdf_reader
#endif
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_io_timer

    contains

        subroutine create_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

#ifdef ENABLE_HDF5
            call create_h5_field_file(basename,     &
                                      overwrite,    &
                                      l_restart)
#endif

#ifdef ENABLE_NETCDF
            call create_netcdf_field_file(basename,     &
                                          overwrite,    &
                                          l_restart)
#endif
        end subroutine create_field_file


        subroutine write_field_step(t, dt)
            double precision,  intent(in) :: t
            double precision,  intent(in) :: dt

            call start_timer(field_io_timer)

#ifdef ENABLE_HDF5
            call write_h5_field_step(t, dt)
#endif

#ifdef ENABLE_NETCDF
            call write_netcdf_field_step(t)
#endif

            call stop_timer(field_io_timer)
        end subroutine write_field_step

        subroutine read_domain(fname, origin, extent, ncells)
            character(*), intent(in)      :: fname
            integer,          intent(out) :: ncells(:)
            double precision, intent(out) :: extent(:), origin(:)

            print *, "read domain"

#ifdef ENABLE_HDF5
            call read_h5_domain(fname, origin, extent, ncells)
#endif

#ifdef ENABLE_NETCDF
!             call read_netcdf_domain(fname, origin, extent, ncells)
#endif
        end subroutine read_domain

end module field_io
