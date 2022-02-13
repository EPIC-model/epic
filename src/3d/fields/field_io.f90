module field_io
    use field_netcdf
    use netcdf_reader
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_io_timer

    contains

        subroutine create_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

            call create_netcdf_field_file(basename,     &
                                          overwrite,    &
                                          l_restart)
        end subroutine create_field_file


        subroutine write_fields(t)
            double precision,  intent(in) :: t

            call start_timer(field_io_timer)

            call write_netcdf_fields(t)
            call stop_timer(field_io_timer)
        end subroutine write_fields

        subroutine read_domain(fname, origin, extent, ncells)
            character(*), intent(in)      :: fname
            integer,          intent(out) :: ncells(:)
            double precision, intent(out) :: extent(:), origin(:)

            call read_netcdf_domain(fname, origin, extent, ncells)
        end subroutine read_domain

end module field_io
