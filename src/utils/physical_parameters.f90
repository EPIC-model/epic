! =============================================================================
! This module contains physical parameters.
! =============================================================================
module physical_parameters
    use constants
    use physical_constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    implicit none

    !FIXME parameters for coriolis and mean wind
!     !non-dimensional ang_vel of earth = t_scale*Omega
!     !t_scale = 142.8571428571 Omega = 7.2921159e-5

    ![m] inverse condensation scale-height
    double precision :: lam_c = 0.001d0

    double precision, protected :: lat_ref

    logical :: l_coriolis = .false.

    ![m/s] angular velocity
    double precision :: ang_vel = twopi / 86400.0d0

    double precision :: lat_degrees = 45.0d0

    ! Coriolis frequency
    double precision, protected :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision, protected :: ft_cor

    ![kg/m**3] saturation specific humidity at ground level
    double precision :: h_0 = 0.015d0

    ![K] mean liquid-water potential temperature
    double precision :: theta_l0 = 300.0d0

    ![] see equation (5) of MPIC paper
    double precision, protected :: glat

    double precision, protected :: glati

    private :: update_physical_constants

    contains

        subroutine update_physical_constants
            glat = gravity * L_v / (c_p * theta_l0)
            glati = one / glat
        end subroutine update_physical_constants

        subroutine read_physical_parameters(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_parameters', grp_ncid)

            ! if no group available --> we use default values
            if (ncerr .ne. 0) then
#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, "WARNING: No physical parameters found! EPIC uses default values."
                endif
#endif
                call update_physical_constants
                return
            endif

            call read_netcdf_attribute_default(grp_ncid, 'specific_humidity', h_0)
            call read_netcdf_attribute_default(grp_ncid, 'liquid_water_potential_temperature', theta_l0)
            call read_netcdf_attribute_default(grp_ncid, 'coriolis', l_coriolis)
            call read_netcdf_attribute_default(grp_ncid, 'angular_velocity', ang_vel)
            call read_netcdf_attribute_default(grp_ncid, 'lat_degrees', lat_degrees)
            call read_netcdf_attribute_default(grp_ncid, 'inverse_condensation_scale_height', lam_c)

            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif

            call update_physical_constants

        end subroutine read_physical_parameters

        subroutine write_physical_parameters(ncid)
            integer, intent(in)     :: ncid
            integer                 :: grp_ncid
            character(*), parameter :: name = 'physical_parameters'

            ncerr = nf90_def_grp(ncid, name, grp_ncid)
            call check_netcdf_error("Faild to create NetCDF group '" // name // "'.")

            call write_netcdf_attribute(grp_ncid, 'specific_humidity', h_0)
            call write_netcdf_attribute(grp_ncid, 'liquid_water_potential_temperature', theta_l0)
            call write_netcdf_attribute(grp_ncid, 'coriolis', l_coriolis)
            call write_netcdf_attribute(grp_ncid, 'angular_velocity', ang_vel)
            call write_netcdf_attribute(grp_ncid, 'lat_degrees', lat_degrees)
            call write_netcdf_attribute(grp_ncid, 'inverse_condensation_scale_height', lam_c)

        end subroutine write_physical_parameters

end module physical_parameters
