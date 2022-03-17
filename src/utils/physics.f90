! =========================================================================================
! This module contains physical quantities.
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
!
!   Dritschel D.G., Böing S.J., Parker D.J., Blyth A.M.
!   The moist parcel-in-cell method for modelling moist convection.
!   Q J R Meteorol Soc. 2018; 144:1695–1718. https://doi.org/10.1002/qj.3319
!
!   Böing S.J., Dritschel D.G., Parker D.J., Blyth, A.M.
!   Comparison of the Moist Parcel-in-Cell (MPIC) model with large-eddy
!   simulation for an idealized cloud.
!   Q J R Meteorol Soc. 2019; 145: 1865– 1881. https://doi.org/10.1002/qj.3532
! =========================================================================================
module physics
    use constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
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

    logical, protected  :: l_coriolis = .false.

    ![°] latitude angle (45° corresponds to standard gravity)
    double precision, protected  :: lat_degrees = 45.0d0

    ![1] MPIC specific: saturation specific humidity at ground level
    double precision, protected  :: q_0 = 0.015d0

    ![m] MPIC specific, scale-height, H
    double precision, protected :: height_c = 1000.0d0

    !
    ! The following quantities are calculated:
    !

    double precision, protected :: lat_ref

    ![1/m] MPIC specific, inverse scale-height
    double precision, protected :: lambda_c

    ! MPIC specific, see equation (5) of MPIC paper
    double precision, protected :: glat

    ! MPIC specific
    double precision, protected :: glati

    ! Coriolis frequency
    double precision, protected :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision, protected :: ft_cor



    ! Warning: This function should only be used by model setup subroutines/programs
    interface set_physical_quantity
        module procedure :: set_physical_quantity_double
        module procedure :: set_physical_quantity_logical
    end interface set_physical_quantity

    interface print_physical_quantity
        module procedure :: print_physical_quantity_double
        module procedure :: print_physical_quantity_integer
        module procedure :: print_physical_quantity_logical
        module procedure :: print_physical_quantity_character
    end interface print_physical_quantity

    private :: update_physical_quantities,          &
               set_physical_quantity_double,        &
               set_physical_quantity_logical,       &
               print_physical_quantity_double,      &
               print_physical_quantity_integer,     &
               print_physical_quantity_logical,     &
               print_physical_quantity_character,   &
               set_coriolis_effects,                &
               set_inverse_scale_height

    contains

        subroutine set_physical_quantity_double(name, val)
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
                case ('latitude_degrees')
                    lat_degrees = val
                case ('saturation_specific_humidity_at_ground_level')
                    q_0 = val
                case ('scale_height')
                    height_c = val
                    call set_inverse_scale_height
                case default
                    print *, "Unknown physical quantity: '" // name // "'."
                    stop
            end select
        end subroutine set_physical_quantity_double

        subroutine set_physical_quantity_logical(name, val)
            character(*), intent(in) :: name
            logical,      intent(in) :: val

            select case (name)
                case ('coriolis')
                    l_coriolis = val
                    call set_coriolis_effects
                case default
                    print *, "Unknown physical quantity: '" // name // "'."
                    stop
            end select
        end subroutine set_physical_quantity_logical

        subroutine set_physical_quantity(name, val)
            character(*),     intent(in) :: name
            double precision, intent(in) :: val

            select case (name)

                case default
                    print *, "Unknown physical constant: '" // name // "'."
                    stop
            end select

        end subroutine set_physical_quantity

        subroutine set_coriolis_effects
            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif
        end subroutine set_coriolis_effects

        subroutine set_inverse_scale_height
            lambda_c = one / height_c
        end subroutine set_inverse_scale_height

        subroutine update_physical_quantities

            glat = gravity * L_v / (c_p * theta_0)
            glati = one / glat

            call set_inverse_scale_height
            call set_coriolis_effects

        end subroutine update_physical_quantities

        subroutine read_physical_quantities(ncid)
            integer, intent(in) :: ncid
            integer             :: grp_ncid

            ncerr = nf90_inq_ncid(ncid, 'physical_quantities', grp_ncid)

            if (ncerr == 0) then
                call read_netcdf_attribute_default(grp_ncid, 'standard_gravity', gravity)
                call read_netcdf_attribute_default(grp_ncid, 'latent_heat_of_vaporization', L_v)
                call read_netcdf_attribute_default(grp_ncid, 'specific_heat', c_p)
                call read_netcdf_attribute_default(grp_ncid, 'temperature_at_sea_level', theta_0)
                call read_netcdf_attribute_default(grp_ncid, 'planetary_angular_velocity', ang_vel)
                call read_netcdf_attribute_default(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
                call read_netcdf_attribute_default(grp_ncid, 'coriolis', l_coriolis)
                call read_netcdf_attribute_default(grp_ncid, 'latitude_degrees', lat_degrees)
                call read_netcdf_attribute_default(grp_ncid, 'scale_height', height_c)
#ifdef ENABLE_VERBOSE
            else
                print *, "WARNING: No physical constants found! EPIC uses default values."
#endif
            endif

            call update_physical_quantities

        end subroutine read_physical_quantities

        subroutine write_physical_quantities(ncid)
            integer, intent(in)     :: ncid
            integer                 :: grp_ncid
            character(*), parameter :: name = 'physical_quantities'

            ncerr = nf90_def_grp(ncid, name, grp_ncid)

            call check_netcdf_error("Faild to create NetCDF group '" // name // "'.")

            call write_netcdf_attribute(grp_ncid, 'standard_gravity', gravity)
            call write_netcdf_attribute(grp_ncid, 'latent_heat_of_vaporization', L_v)
            call write_netcdf_attribute(grp_ncid, 'specific_heat', c_p)
            call write_netcdf_attribute(grp_ncid, 'temperature_at_sea_level', theta_0)
            call write_netcdf_attribute(grp_ncid, 'planetary_angular_velocity', ang_vel)
            call write_netcdf_attribute(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
            call write_netcdf_attribute(grp_ncid, 'coriolis', l_coriolis)
            call write_netcdf_attribute(grp_ncid, 'coriolis_frequency', f_cor)
            call write_netcdf_attribute(grp_ncid, 'planetary_vorticity', ft_cor)
            call write_netcdf_attribute(grp_ncid, 'latitude_degrees', lat_degrees)
            call write_netcdf_attribute(grp_ncid, 'scale_height', height_c)
            call write_netcdf_attribute(grp_ncid, 'inverse_scale_height', lambda_c)

        end subroutine write_physical_quantities

        subroutine print_physical_quantities
            write(*, "(a)") 'List of physical quantities (in MKS units):'
            write(*, "(a)") repeat("-", 78)
            call print_physical_quantity('standard gravity', gravity, 'm/s^2')
            call print_physical_quantity('latent heat of vaporization', L_v, 'J/kg')
            call print_physical_quantity('specific heat', c_p, 'J/(kg*K)')
            call print_physical_quantity('temperature at sea level', theta_0, 'K')
            call print_physical_quantity('planetary angular velocity', ang_vel, 'rad/s')
            call print_physical_quantity('saturation specific humidity at ground level', q_0)
            call print_physical_quantity('coriolis', l_coriolis)
            call print_physical_quantity('coriolis frequency', f_cor, '1/s')
            call print_physical_quantity('planetary vorticity', ft_cor, '1/s')
            call print_physical_quantity('latitude degrees', lat_degrees, 'deg')
            call print_physical_quantity('scale height', height_c, 'm')
            call print_physical_quantity('inverse scale height', lambda_c, '1/m')
            write(*, *) ''
        end subroutine print_physical_quantities

        subroutine print_physical_quantity_double(name, val, unit)
            character(*),           intent(in) :: name
            double precision,       intent(in) :: val
            character(*), optional, intent(in) :: unit
            character(64)                      :: fix_length_name

            fix_length_name = name
            if (present(unit)) then
                fix_length_name = name // ', ' // unit
            endif
            write (*, "(a, 1p,e14.7)") fix_length_name, val
        end subroutine print_physical_quantity_double

        subroutine print_physical_quantity_integer(name, val, unit)
            character(*),           intent(in) :: name
            integer,                intent(in) :: val
            character(*), optional, intent(in) :: unit
            character(64)                      :: fix_length_name

            fix_length_name = name
            if (present(unit)) then
                fix_length_name = name // ', ' // unit
            endif
            write (*, "(a, I14)") fix_length_name, val
        end subroutine print_physical_quantity_integer

        subroutine print_physical_quantity_logical(name, val, unit)
            character(*),           intent(in) :: name
            logical,                intent(in) :: val
            character(*), optional, intent(in) :: unit

            if (val) then
                call print_physical_quantity_character(name, 'true')
            else
                call print_physical_quantity_character(name, 'false')
            endif
        end subroutine print_physical_quantity_logical

        subroutine print_physical_quantity_character(name, val, unit)
            character(*),           intent(in) :: name
            character(*),           intent(in) :: val
            character(*), optional, intent(in) :: unit
            character(64)                      :: fix_length_name

            fix_length_name = name
            write (*, "(a, a14)") fix_length_name, val
        end subroutine print_physical_quantity_character

end module physics
