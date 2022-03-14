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
    double precision, protected :: lam_c = 0.001d0

    ![] see equation (5) of MPIC paper
    double precision, protected :: glat

    double precision, protected :: glati

    double precision, protected :: lat_ref

    logical, protected :: l_coriolis

    ![m/s] angular velocity
    double precision, protected :: ang_vel

    double precision, protected :: lat_degrees

    ! Coriolis frequency
    double precision, protected :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision, protected :: ft_cor

    contains

        subroutine update_physical_parameters

             glat = gravity * L_v / (c_p * theta_l0)
             glati = one / glat

        end subroutine update_physical_parameters

        subroutine read_physical_parameters(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_parameters', grp_ncid)

            ! if no group available, leave the function
            if (ncerr .ne. 0) then
#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, "WARNING: No physical parameters found! EPIC uses default values."
                endif
#endif
                return
            endif

            call read_netcdf_attribute_default(grp_ncid, 'coriolis', l_coriolis, .false.)
            call read_netcdf_attribute_default(grp_ncid, 'angular_velocity', ang_vel, twopi / 86400.0d0)
            call read_netcdf_attribute_default(grp_ncid, 'lat_degrees', lat_degrees, 45.0d0)


            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif

        end subroutine read_physical_parameters

        subroutine write_physical_parameters

        end subroutine write_physical_parameters

end module physical_parameters
