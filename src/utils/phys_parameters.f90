! =============================================================================
! This module contains modifiable physical parameters.
! =============================================================================
module phys_parameters
    use constants
!     use options, only : l_coriolis, lat_degrees, ang_vel
    use phys_constants
    implicit none

    ![m] inverse condensation scale-height
    double precision, parameter :: lam_c = 0.001d0

    ![] see equation (5) of MPIC paper
    double precision, parameter :: glat = gravity * L_v / (c_p * theta_l0)

    double precision, parameter :: glati = one / glat

    !FIXME comment
    double precision :: lat_ref

    ! Coriolis frequency
    double precision :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision :: ft_cor

    private :: get_phys_parameter

    contains

        subroutine update_phys_parameters
!             if (l_coriolis) then
!                 lat_ref = lat_degrees * deg2rad
!                 f_cor  = two * ang_vel * dsin(lat_ref)
!                 ft_cor = two * ang_vel * dcos(lat_ref)
!             else
!                 f_cor  = zero
!                 ft_cor = zero
!             endif
        end subroutine update_phys_parameters

        subroutine get_phys_parameter(ncid, name, val)
            integer,      intent(in)        :: ncid
            character(*), intent(in)        :: name
            double precision, intent(inout) :: val      ! needs to be "inout" if parameter not available

            if (has_attribute(grp_ncid, name)) then
                call read_netcdf_attribute(grp_ncid, name, val)
#ifdef ENABLE_VERBOSE
                    if (verbose) then
                        print *, "Found physical input parameter '" // name // "'."
                    endif
            else
                if (verbose) then
                    print *, "WARNING: Using default value of '" // name // "'."
                endif
#endif
            endif

        end subroutine get_phys_parameter

        subroutine read_phys_parameters(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_parameters', grp_ncid)

            ! if no group available, leave the function
            if (ncerr .ne. 0) then
#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, "No physical input parameters found! EPIC uses default values."
                endif
#endif
                return
            endif

            call get_phys_parameter(grp_ncid, 'gravity', gravity)
            call get_phys_parameter(grp_ncid, 'latent_heat', L_v)
            call get_phys_parameter(grp_ncid, 'specific_heat', c_p)
            call get_phys_parameter(grp_ncid, 'specific_humidity', h_0)
            call get_phys_parameter(grp_ncid, 'liquid_water_potential_temperature', theta_l0)

        end subroutine read_phys_parameters

end module phys_parameters
