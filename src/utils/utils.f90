module utils
    use constants, only : one, ndim
    use field_hdf5
    use parcel_hdf5
    use parcel_diagnostics
    use parcel_container, only : n_parcels
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
#ifndef NDEBUG
    use parcel_interpl, only : vol2grid_symmetry_error
#endif
    use field_diagnostics, only : write_h5_field_stats_step

    implicit none

    integer :: nfw  = 0    ! number of field writes to h5
    integer :: npw  = 0    ! number of parcel writes to h5
    integer :: nspw = 0    ! number of parcel diagnostics writes to h5
    integer :: nsfw = 0    ! number of field diagnostics writes to h5

    private :: nfw, npw, nspw, nsfw

    contains

        ! Write last step to the H5 files. For the time step dt, it
        ! writes zero.
        ! @param[in] t is the time
        subroutine write_last_step(t)
            double precision,  intent(in) :: t
            double precision              :: velocity(n_parcels, ndim)
            double precision              :: strain(n_parcels, 4) !FIXME
            double precision              :: vorticity(n_parcels)

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(tbuoyg, vtend)

            call grid2par(velocity, vorticity, strain)

            call calc_parcel_diagnostics(velocity)

            call write_step(t, zero, .true.)
        end subroutine write_last_step

        ! Write step to the H5 files.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        ! @param[in] l_force a logical to force a write (optional)
        subroutine write_step(t, dt, l_force)
            use options, only : output
            double precision,  intent(in) :: t
            double precision,  intent(in) :: dt
            logical, optional, intent(in) :: l_force
            double precision              :: neg = one
#ifndef NDEBUG
            logical                      :: do_vol2grid_sym_err = .true.
#endif

            if (present(l_force)) then
                if (l_force) then
                    neg = -one
                endif
            endif

            ! make sure we always write initial setup
            if (output%h5_write_fields .and. &
                (t + epsilon(zero) >= neg * dble(nfw) * output%h5_field_freq)) then
#ifndef NDEBUG
                call vol2grid_symmetry_error
                do_vol2grid_sym_err = .false.
#endif
                call write_h5_field_step(nfw, t, dt)
            endif


            if (output%h5_write_parcels .and. &
                (t + epsilon(zero) >= neg * dble(npw) * output%h5_parcel_freq)) then
                call write_h5_parcel_step(npw, t, dt)
            endif

            if (output%h5_write_parcel_stats .and. &
                (t + epsilon(zero) >= neg * dble(nspw) * output%h5_parcel_stats_freq)) then
                call write_h5_parcel_stats_step(nspw, t, dt)
            endif

            if (output%h5_write_field_stats .and. &
                (t + epsilon(zero) >= neg * dble(nsfw) * output%h5_field_stats_freq)) then

#ifndef NDEBUG
                if (do_vol2grid_sym_err) then
                    call vol2grid_symmetry_error
                endif
#endif
                call write_h5_field_stats_step(nsfw, t, dt)
            endif
        end subroutine write_step

end module utils
