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
    implicit none

    ![m/s**2]
    double precision, protected :: gravity

    ![J/kg] latent heat of condensation
    double precision, protected :: L_v

    ![J/(kg*K)] specific heat at constant pressure
    double precision, protected :: c_p

    ![kg/m**3] saturation specific humidity at ground level
    double precision, protected :: h_0

    ![K] mean liquid-water potential temperature
    double precision, protected :: theta_l0

    contains

        subroutine read_physical_constants(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_constants', grp_ncid)

            ! if no group available, leave the function
            if (ncerr .ne. 0) then
#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, "WARNING: No physical constants found! EPIC uses default values."
                endif
#endif
                return
            endif

            call read_netcdf_attribute_default(grp_ncid, 'gravity', gravity, 9.81d0)
            call read_netcdf_attribute_default(grp_ncid, 'latent_heat', L_v, 2.501e6)
            call read_netcdf_attribute_default(grp_ncid, 'specific_heat', c_p, 1005.0d0)
            call read_netcdf_attribute_default(grp_ncid, 'specific_humidity', h_0, 0.015d0)
            call read_netcdf_attribute_default(grp_ncid, 'liquid_water_potential_temperature', theta_l0, 300.0d0)

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
            call write_netcdf_attribute(grp_ncid, 'specific_humidity', h_0)
            call write_netcdf_attribute(grp_ncid, 'liquid_water_potential_temperature', theta_l0)

        end subroutine write_physical_constants

end module physical_constants
