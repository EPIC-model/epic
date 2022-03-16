! =========================================================================================
! This module contains physical constants.
!
! Note: We cannot use the keyword "parameter" since we read in the values
!       of the quantities from the NetCDF file. Hence, they are not known
!       at compile-time which is required for "parameter variables". However,
!       the keyword "protected" (Fortran2003) allows to make the quantities
!       non-modifiable outside the module (14 March 2022,
!       https://stackoverflow.com/questions/15020460/protected-global-variables-in-fortran)
!
! References:
!   American Meteorological Society, 2022: Standard gravity. Glossary of Meteorology,
!   https://glossary.ametsoc.org/wiki/Standard_gravity
!
!   American Meteorological Society, 2022: Latent heat. Glossary of Meteorology,
!   https://glossary.ametsoc.org/wiki/Latent_heat
!
!   American Meteorological Society, 2022: Specific heat capacity. Glossary of Meteorology,
!   https://glossary.ametsoc.org/wiki/Specific_heat_capacity
!
!   American Meteorological Society, 2022: Standard atmosphere. Glossary of Meteorology,
!   https://glossary.ametsoc.org/wiki/Standard_atmosphere
!
!   American Meteorological Society, 2022: Angular velocity of the earth. Glossary of Meteorology,
!   https://glossary.ametsoc.org/wiki/Angular_velocity_of_the_earth
! =========================================================================================
module physical_constants
    use constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use iomanip, only : print_quantity
    implicit none

    ![m/s**2] standard gravity (i.e. at 45° latitude and mean sea level):
    double precision, protected :: gravity = 9.80616d0

    ![J/kg] latent heat of vaporization:
    double precision, protected :: L_v = 2.501e6

    ![J/(kg*K)] specific heat at constant pressure and dry air:
    double precision, protected :: c_p = 1005.7d0

    ![K] temperature at zero pressure altitude (15°C)
    double precision, protected :: theta_0 = 288.15d0

    ![rad/s] default value: angular velocity of the earth
    double precision :: ang_vel = 0.000072921d0

    contains

        ! Warning: This function should only be used by model setup subroutines/programs
        subroutine set_physical_constant(name, val)
            character(*),     intent(in) :: name
            double precision, intent(in) :: val

            select case (name)
                case ('standard_gravity')
                    gravity = val
                case ('latent_heat_of_vaporization')
                    L_v = val
                case ('specific_heat')
                    c_p = val
                case ('temperature_at_sea_level')
                    theta_0 = val
                case ('planetary_angular_velocity')
                    ang_vel = val
                case default
                    print *, "Unknown physical constant: '", trim(name), "'."
                    stop
            end select

        end subroutine set_physical_constant

        subroutine read_physical_constants(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_constants', grp_ncid)

            if (ncerr == 0) then
                call read_netcdf_attribute_default(grp_ncid, 'standard_gravity', gravity)
                call read_netcdf_attribute_default(grp_ncid, 'latent_heat_of_vaporization', L_v)
                call read_netcdf_attribute_default(grp_ncid, 'specific_heat', c_p)
                call read_netcdf_attribute_default(grp_ncid, 'temperature_at_sea_level', theta_0)
                call read_netcdf_attribute_default(grp_ncid, 'planetary_angular_velocity', ang_vel)
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

            call write_netcdf_attribute(grp_ncid, 'standard_gravity', gravity)
            call write_netcdf_attribute(grp_ncid, 'latent_heat_of_vaporization', L_v)
            call write_netcdf_attribute(grp_ncid, 'specific_heat', c_p)
            call write_netcdf_attribute(grp_ncid, 'temperature_at_sea_level', theta_0)
            call write_netcdf_attribute(grp_ncid, 'planetary_angular_velocity', ang_vel)

        end subroutine write_physical_constants

        subroutine print_physical_constants
            write(*, "(a)") 'List of physical constants (in MKS units):'
            write(*, "(a)") repeat("-", 78)
            call print_quantity('standard gravity', gravity, 'm/s^2')
            call print_quantity('latent heat of vaporization', L_v, 'J/kg')
            call print_quantity('specific heat', c_p, 'J/(kg*K)')
            call print_quantity('temperature at sea level', theta_0, 'K')
            call print_quantity('planetary angular velocity', pl_ang_vel, 'rad/s')
            write(*, *) ''
        end subroutine print_physical_constants

end module physical_constants
