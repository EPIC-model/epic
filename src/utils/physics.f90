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
    use iso_fortran_env, only : IOSTAT_END
    use ape_density, only : l_ape_density
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

    logical, protected  :: l_planet_vorticity = .false.

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

    ! planetary vorticity (all three components)
    double precision, protected :: f_cor(3)

    ! domain-averaged potential energy reference
    double precision :: peref

    ! is .true. when 'peref' was read in
    logical, protected :: l_peref

    ! 'none', 'sorting' or 'ape density'
    character(len=11) :: ape_calculation = "sorting"

    interface print_physical_quantity
        module procedure :: print_physical_quantity_double
        module procedure :: print_physical_quantity_integer
        module procedure :: print_physical_quantity_logical
        module procedure :: print_physical_quantity_character
    end interface print_physical_quantity

    private :: update_physical_quantities,          &
               print_physical_quantity_double,      &
               print_physical_quantity_integer,     &
               print_physical_quantity_logical,     &
               print_physical_quantity_character

    contains

        ! This subroutine is only used in model setups
        subroutine read_physical_quantities_from_namelist(fname)
            character(*), intent(in) :: fname
            integer                  :: ios
            integer                  :: fn = 2
            logical                  :: exists = .false.

            ! namelist definitions
            namelist /PHYSICS/ gravity,             &
                               L_v,                 &
                               c_p,                 &
                               theta_0,             &
                               ang_vel,             &
                               l_planet_vorticity,  &
                               lat_degrees,         &
                               q_0,                 &
                               height_c

            ! check whether file exists
            inquire(file=fname, exist=exists)

            if (exists .eqv. .false.) then
                print *, 'Error: input file "' // fname // '" does not exist.'
                stop
            endif

            ! open and read Namelist file.
            open(action='read', file=fname, iostat=ios, newunit=fn)

            read(nml=PHYSICS, iostat=ios, unit=fn)

            if (ios == IOSTAT_END) then
                ! physical constants/parameters not present
            else if (ios /= 0) then
                print *, 'Error: invalid Namelist format.'
                stop
            endif

            close(fn)

            call update_physical_quantities

        end subroutine read_physical_quantities_from_namelist

        subroutine update_physical_quantities

            glat = gravity * L_v / (c_p * theta_0)
            glati = one / glat

            lambda_c = one / height_c

            if (l_planet_vorticity) then
                lat_ref = lat_degrees * deg2rad
                f_cor(1) = zero
                f_cor(2) = two * ang_vel * dcos(lat_ref)
                f_cor(3) = two * ang_vel * dsin(lat_ref)
            else
                f_cor = zero
            endif
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
                call read_netcdf_attribute_default(grp_ncid, 'planetary_vorticity', l_planet_vorticity)
                call read_netcdf_attribute_default(grp_ncid, 'latitude_degrees', lat_degrees)
                call read_netcdf_attribute_default(grp_ncid, 'scale_height', height_c)
                call read_netcdf_attribute_default(grp_ncid, 'ape_calculation', ape_calculation)

                l_peref = .false.
                select case (trim(ape_calculation))
                    case ('sorting')
                        l_peref = has_attribute(grp_ncid, 'reference_potential_energy')
                        if (l_peref) then
                            print *, "Found float attribute 'reference_potential_energy'."
                            call read_netcdf_attribute(grp_ncid, 'reference_potential_energy', peref)
                        else
                            print *, "No float attribute 'reference_potential_energy'. It will be computed."
                        endif
                    case ('ape density')
                        if (.not. l_ape_density) then
                            print *, "In order to use the APE calculation, you must provide"
                            print *, "the APE density function in utils/ape_density.f90"
                            error stop
                        endif
                            print *, "APE calculation using APE density function."
                    case default
                        print *, "WARNING: No APE calculated!"
                        ape_calculation = 'none'
                    end select


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

            call check_netcdf_error("Failed to create NetCDF group '" // name // "'.")

            call write_netcdf_attribute(grp_ncid, 'standard_gravity', gravity)
            call write_netcdf_attribute(grp_ncid, 'latent_heat_of_vaporization', L_v)
            call write_netcdf_attribute(grp_ncid, 'specific_heat', c_p)
            call write_netcdf_attribute(grp_ncid, 'temperature_at_sea_level', theta_0)
            call write_netcdf_attribute(grp_ncid, 'planetary_angular_velocity', ang_vel)
            call write_netcdf_attribute(grp_ncid, 'saturation_specific_humidity_at_ground_level', q_0)
            call write_netcdf_attribute(grp_ncid, 'planet_vorticity', l_planet_vorticity)
            call write_netcdf_attribute(grp_ncid, 'latitude_degrees', lat_degrees)
            call write_netcdf_attribute(grp_ncid, 'scale_height', height_c)
            if (l_peref) then
                call write_netcdf_attribute(grp_ncid, 'reference_potential_energy', peref)
            endif
            call write_netcdf_attribute(grp_ncid, 'ape_calculation', ape_calculation)

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
            call print_physical_quantity('planetary vorticity', l_planet_vorticity)
            call print_physical_quantity('vertical planetary vorticity', f_cor(3), '1/s')
            call print_physical_quantity('horizontal planetary vorticity', f_cor(2), '1/s')
            call print_physical_quantity('latitude degrees', lat_degrees, 'deg')
            call print_physical_quantity('scale height', height_c, 'm')
            call print_physical_quantity('inverse scale height', lambda_c, '1/m')
            if (l_peref) then
                call print_physical_quantity('reference potential energy', peref, 'm^2/s^2')
            endif
            call print_physical_quantity('APE calculation', ape_calculation)
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
                call print_physical_quantity_character(name, 'true', unit)
            else
                call print_physical_quantity_character(name, 'false', unit)
            endif
        end subroutine print_physical_quantity_logical

        subroutine print_physical_quantity_character(name, val, unit)
            character(*),           intent(in) :: name
            character(*),           intent(in) :: val
            character(*), optional, intent(in) :: unit
            character(64)                      :: fix_length_name

            fix_length_name = name
            if (present(unit)) then
                fix_length_name = name // ', ' // unit
            endif
            write (*, "(a, a14)") fix_length_name, val
        end subroutine print_physical_quantity_character

end module physics
