module field_io
#ifdef ENABLE_NETCDF
    use field_netcdf
    use netcdf_reader
#else
    use field_hdf5
    use h5_reader
#endif
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_io_timer

    contains

        subroutine create_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

#ifdef ENABLE_NETCDF
            call create_netcdf_field_file(basename,     &
                                          overwrite,    &
                                          l_restart)
#else
            call create_h5_field_file(basename,     &
                                      overwrite,    &
                                      l_restart)
#endif
        end subroutine create_field_file


        subroutine write_fields(t)
            double precision,  intent(in) :: t

            call start_timer(field_io_timer)

#ifdef ENABLE_NETCDF
            call write_netcdf_fields(t)
#else
            call write_h5_fields(t)
#endif
            call stop_timer(field_io_timer)
        end subroutine write_fields

        subroutine read_domain(fname, origin, extent, ncells)
            character(*), intent(in)      :: fname
            integer,          intent(out) :: ncells(:)
            double precision, intent(out) :: extent(:), origin(:)

#ifdef ENABLE_NETCDF
            call read_netcdf_domain(fname, origin, extent, ncells)
#else
            call read_h5_domain(fname, origin, extent, ncells)
#endif
        end subroutine read_domain

end module field_io
