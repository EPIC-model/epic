module diagnostics
    use parcel_container
    use parameters, only : vcell
    use options, only : grid
    use fields
    use hdf5
    use writer, only : h5file,                        &
                       h5err,                         &
                       write_h5_double_scalar_attrib, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none


    contains

        function get_rms_volume_error() result(rms)
            double precision              :: n
            double precision              :: rms ! rms volume error
            double precision, allocatable :: V(:, :, :)

            ! remove halo cells
            V = volume_f(1:grid(1), 1:grid(2), :)

            n = size(V)
            rms = sqrt(sum((V - vcell) ** 2) / n) / vcell
        end function get_rms_volume_error


        subroutine write_h5_diagnostics(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            integer(hid_t)                :: step_group
            character(:), allocatable     :: step
            character(:), allocatable     :: name
            double precision              :: rms_v

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/diagnostics"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            !
            ! write diagnostics
            !

            rms_v = get_rms_volume_error()
            call write_h5_double_scalar_attrib(group, "rms volume error", rms_v)


            ! close all
            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)
        end subroutine write_h5_diagnostics

end module diagnostics
