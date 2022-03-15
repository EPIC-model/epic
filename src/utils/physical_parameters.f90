! =============================================================================
! This module contains physical parameters.
! =============================================================================
module physical_parameters
    use constants
    use physical_constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use iomanip, only : print_key_value_pair
    implicit none

    !FIXME parameters for coriolis and mean wind
!     !non-dimensional ang_vel of earth = t_scale*Omega
!     !t_scale = 142.8571428571 Omega = 7.2921159e-5

    logical :: l_coriolis = .false.

    ![m/s] angular velocity
    double precision :: ang_vel = twopi / 86400.0d0

    double precision :: lat_degrees = 45.0d0

    ![kg/m**3] saturation specific humidity at ground level
    double precision :: q_0 = 0.015d0

    ![K] reference virtual/density temperature
    double precision :: theta_v0 = 298.6268656716418d0

    ![] see equation (5) of MPIC paper
    double precision, protected :: glat

    double precision, protected :: glati

    double precision, protected :: lat_ref

    ! Coriolis frequency
    double precision, protected :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision, protected :: ft_cor

    ![m] inverse condensation scale-height
    double precision :: lam_c

    private :: update_physical_parameters

    contains

        subroutine update_physical_parameters
            glat = gravity * L_v / (c_p * theta_v0)
            glati = one / glat

            lam_c = 1000.0d0 * (10.0d0 / gravity)
            lam_c = one / lam_c
        end subroutine update_physical_parameters

        subroutine read_physical_parameters(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_parameters', grp_ncid)

            if (ncerr == 0) then
                call read_netcdf_attribute_default(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
                call read_netcdf_attribute_default(grp_ncid, 'reference_virtual_temperature', theta_v0)
                call read_netcdf_attribute_default(grp_ncid, 'coriolis', l_coriolis)
                call read_netcdf_attribute_default(grp_ncid, 'angular_velocity', ang_vel)
                call read_netcdf_attribute_default(grp_ncid, 'latitude_degrees', lat_degrees)
#ifdef ENABLE_VERBOSE
            else
                print *, "WARNING: No physical parameters found! EPIC uses default values."
#endif
            endif


            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif

            call update_physical_parameters

        end subroutine read_physical_parameters

        subroutine write_physical_parameters(ncid)
            integer, intent(in)     :: ncid
            integer                 :: grp_ncid
            character(*), parameter :: name = 'physical_parameters'

            ncerr = nf90_def_grp(ncid, name, grp_ncid)
            call check_netcdf_error("Faild to create NetCDF group '" // name // "'.")

            call write_netcdf_attribute(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
            call write_netcdf_attribute(grp_ncid, 'reference_virtual_temperature', theta_v0)
            call write_netcdf_attribute(grp_ncid, 'coriolis', l_coriolis)
            call write_netcdf_attribute(grp_ncid, 'coriolis_frequency', f_cor)
            call write_netcdf_attribute(grp_ncid, 'planetary_vorticity', ft_cor)
            call write_netcdf_attribute(grp_ncid, 'angular_velocity', ang_vel)
            call write_netcdf_attribute(grp_ncid, 'latitude_degrees', lat_degrees)
            call write_netcdf_attribute(grp_ncid, 'inverse_condensation_scale_height', lam_c)

        end subroutine write_physical_parameters

        subroutine print_physical_parameters
            write(*, "(a)") 'List of physical parameters (in MKS units):'
            write(*, "(a)") repeat("-", 78)
            call print_key_value_pair('saturation specific humidity at ground level', q_0)
            call print_key_value_pair('reference virtual temperature', theta_v0)
            call print_key_value_pair('coriolis', l_coriolis)
            call print_key_value_pair('coriolis frequency', f_cor)
            call print_key_value_pair('planetary vorticity', ft_cor)
            call print_key_value_pair('angular velocity', ang_vel)
            call print_key_value_pair('latitude degrees', lat_degrees)
            call print_key_value_pair('inverse condensation scale height', lam_c)
            write(*, *) ''
        end subroutine print_physical_parameters

end module physical_parameters
