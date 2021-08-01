! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower
    use parcel_container, only : parcels, n_parcels
    use h5_utils
    use h5_writer
    use omp_lib
    use timer, only : start_timer, stop_timer
    implicit none


    private

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    ! peref : potential energy reference
    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: peref, pe, ke

    integer :: hdf5_parcel_stat_timer

    public :: create_h5_parcel_stat_file,   &
              init_parcel_diagnostics,     &
              write_h5_parcel_stats_step,  &
              hdf5_parcel_stat_timer

    contains

        subroutine create_h5_parcel_stat_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical,      intent(in) :: overwrite

            h5fname =  basename // '_parcel_diagnostics.hdf5'

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_char_scalar_attrib(h5file_id, 'output_type', 'parcel diagnostics')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id)

            call close_h5_file(h5file_id)

        end subroutine create_h5_parcel_stat_file

        ! compute the reference potential energy
        subroutine init_parcel_diagnostics
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = one / extent(1)
            zmean = lower(2) + f12 * gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * parcels%volume(ii(n))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo
        end subroutine init_parcel_diagnostics


        subroutine calculate_diagnostics
            integer          :: n
            double precision :: b, z, vel(2), vol, zmin

            ! reset
            ke = zero
            pe = zero

            zmin = lower(2)

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z) reduction(+: ke, pe)
            do n = 1, n_parcels

                vel = parcels%velocity(n, :)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(n, 2) - zmin

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                ! potential energy
                pe = pe - b * z * vol
            enddo
            !$omp end do
            !$omp end parallel

            ke = f12 * ke
            pe = pe - peref
        end subroutine calculate_diagnostics


        subroutine write_h5_parcel_stats_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            integer(hid_t)                  :: group
            character(:), allocatable       :: name
            logical                         :: created

            call start_timer(hdf5_parcel_stat_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write parcel diagnostics to h5"
            endif
#endif

            call calculate_diagnostics

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif


            call write_h5_double_scalar_attrib(group, "t", t)

            call write_h5_double_scalar_attrib(group, "dt", dt)

            !
            ! write diagnostics
            !
            call write_h5_double_scalar_attrib(group, "potential energy", pe)
            call write_h5_double_scalar_attrib(group, "kinetic energy", ke)
            call write_h5_double_scalar_attrib(group, "total energy", ke + pe)
            call write_h5_int_scalar_attrib(group, "num parcel", n_parcels)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_group(group)

            call close_h5_file(h5file_id)

            call stop_timer(hdf5_parcel_stat_timer)

        end subroutine write_h5_parcel_stats_step
end module parcel_diagnostics
