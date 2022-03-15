! =========================================================================================
! This module contains physical constants.
!
! Note: We cannot use the keyword "parameter" since we read in the values
!       of the quantities from the NetCDF file. Hence, they are not known
!       at compile-time which is required for "parameter variables". However,
!       the keyword "protected" (Fortran2003) allows to make the quantities
!       non-modifiable outside the module (14 March 2022,
!       https://stackoverflow.com/questions/15020460/protected-global-variables-in-fortran)
! =========================================================================================
module physical_constants
    use constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use iomanip, only : print_key_value_pair
    implicit none

    ![m/s**2]
    double precision, protected :: gravity = 9.81d0

    ![J/kg] latent heat of condensation
    double precision, protected :: L_v = 2.501e6

    ![J/(kg*K)] specific heat at constant pressure
    double precision, protected :: c_p = 1005.0d0

    contains

        subroutine read_physical_constants(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_constants', grp_ncid)

            if (ncerr == 0) then
                call read_netcdf_attribute_default(grp_ncid, 'gravity', gravity)
                call read_netcdf_attribute_default(grp_ncid, 'latent_heat', L_v)
                call read_netcdf_attribute_default(grp_ncid, 'specific_heat', c_p)
#ifdef ENABLE_VERBOSE
            else
                print *, "WARNING: No physical constants found! EPIC uses default values."
#endif
            endif

        end subroutine read_physical_constants

        subroutine write_physical_constants(ncid)
            integer, intent(in)     :: ncid
            integer                 :: grp_ncid
            character(*), parameter :: name = 'physical_constants'

            ncerr = nf90_def_grp(ncid, name, grp_ncid)

            call check_netcdf_error("Faild to create NetCDF group '" // name // "'.")

            call write_netcdf_attribute(grp_ncid, 'gravity', gravity)
            call write_netcdf_attribute(grp_ncid, 'latent_heat', L_v)
            call write_netcdf_attribute(grp_ncid, 'specific_heat', c_p)

        end subroutine write_physical_constants

        subroutine print_physical_constants
            write(*, "(a)") 'List of physical constants:'
            write(*, "(a)") repeat("-", 78)
            call print_key_value_pair('gravity', gravity)
            call print_key_value_pair('latent heat', L_v)
            call print_key_value_pair('specific heat', c_p)
            write(*, *) ''
        end subroutine print_physical_constants

end module physical_constants
