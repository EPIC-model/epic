! =============================================================================
! This module contains physical parameters.
!
! References:
!   Dritschel D.G., Böing S.J., Parker D.J., Blyth A.M.
!   The moist parcel-in-cell method for modelling moist convection.
!   Q J R Meteorol Soc. 2018; 144:1695–1718. https://doi.org/10.1002/qj.3319
!
!   Böing S.J., Dritschel D.G., Parker D.J., Blyth, A.M.
!   Comparison of the Moist Parcel-in-Cell (MPIC) model with large-eddy
!   simulation for an idealized cloud.
!   Q J R Meteorol Soc. 2019; 145: 1865– 1881. https://doi.org/10.1002/qj.3532
! =============================================================================
module physical_parameters
    use constants
    use physical_constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use iomanip, only : print_quantity
    implicit none

    logical, protected  :: l_coriolis = .false.

    ![°] latitude angle (45° corresponds to standard gravity)
    double precision, protected  :: lat_degrees = 45.0d0

    ![1] MPIC specific: saturation specific humidity at ground level
    double precision, protected  :: q_0 = 0.015d0

    ! MPIC specific, see equation (5) of MPIC paper
    double precision, protected :: glat

    ! MPIC specific
    double precision, protected :: glati

    double precision, protected :: lat_ref

    ! Coriolis frequency
    double precision, protected :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision, protected :: ft_cor

    ![m] MPIC specific, scale-height, H
    double precision, protected :: height_c = 1000.0d0

    ![1/m] MPIC specific, inverse scale-height
    double precision, protected :: lambda_c

    interface set_physical_parameter
        module procedure :: set_physical_parameter_double
        module procedure :: set_physical_parameter_logical
    end interface set_physical_parameter

    private :: update_physical_parameters,      &
               set_physical_parameter_double,   &
               set_physical_parameter_logical

    contains

        subroutine set_physical_parameter_double(name, val)
            character(*),     intent(in) :: name
            double precision, intent(in) :: val

            select case (name)
                case ('latitude_degrees')
                    lat_degrees = val
                case ('saturation_specific_humidity_at_ground_level')
                    q_0 = val
                case ('scale_height')
                    height_c = val
                case default
                    print *, "Unknown physical parameter: '" // name // "'."
                    stop
            end select
        end subroutine set_physical_parameter_double

        subroutine set_physical_parameter_logical(name, val)
            character(*), intent(in) :: name
            logical,      intent(in) :: val

            select case (name)
                case ('coriolis')
                    l_coriolis = val
                case default
                    print *, "Unknown physical parameter: '" // name // "'."
                    stop
            end select
        end subroutine set_physical_parameter_logical

        subroutine update_physical_parameters
            glat = gravity * L_v / (c_p * theta_0)
            glati = one / glat

            lambda_c = one / height_c

            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif
        end subroutine update_physical_parameters

        subroutine read_physical_parameters(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_parameters', grp_ncid)

            if (ncerr == 0) then
                call read_netcdf_attribute_default(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
                call read_netcdf_attribute_default(grp_ncid, 'coriolis', l_coriolis)
                call read_netcdf_attribute_default(grp_ncid, 'latitude_degrees', lat_degrees)
                call read_netcdf_attribute_default(grp_ncid, 'scale_height', height_c)
#ifdef ENABLE_VERBOSE
            else
                print *, "WARNING: No physical parameters found! EPIC uses default values."
#endif
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
            call write_netcdf_attribute(grp_ncid, 'coriolis', l_coriolis)
            call write_netcdf_attribute(grp_ncid, 'coriolis_frequency', f_cor)
            call write_netcdf_attribute(grp_ncid, 'planetary_vorticity', ft_cor)
            call write_netcdf_attribute(grp_ncid, 'latitude_degrees', lat_degrees)
            call write_netcdf_attribute(grp_ncid, 'scale_height', height_c)
            call write_netcdf_attribute(grp_ncid, 'inverse_scale_height', lambda_c)

        end subroutine write_physical_parameters

        subroutine print_physical_parameters
            write(*, "(a)") 'List of physical parameters (in MKS units):'
            write(*, "(a)") repeat("-", 78)
            call print_quantity('saturation specific humidity at ground level', q_0)
            call print_quantity('coriolis', l_coriolis)
            call print_quantity('coriolis frequency', f_cor, '1/s')
            call print_quantity('planetary vorticity', ft_cor, '1/s')
            call print_quantity('latitude degrees', lat_degrees, 'deg')
            call print_quantity('scale height', height_c, 'm')
            call print_quantity('inverse scale height', lambda_c, '1/m')
            write(*, *) ''
        end subroutine print_physical_parameters

end module physical_parameters
