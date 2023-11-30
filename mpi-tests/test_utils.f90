module test_utils
    use mpi_timer
    use parcel_container, only : resize_timer
    use parcel_split_mod, only : split_timer
    use parcel_merging, only : merge_timer
    use parcel_nearest, only : merge_nearest_timer, merge_tree_resolve_timer
    use parcel_correction, only : lapl_corr_timer,        &
                                  grad_corr_timer,        &
                                  vort_corr_timer
    use parcel_diagnostics, only : parcel_stats_timer
    use parcel_netcdf, only : parcel_io_timer
    use parcel_diagnostics_netcdf, only : parcel_stats_io_timer
    use field_netcdf, only : field_io_timer
    use field_diagnostics, only : field_stats_timer
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use inversion_mod, only : vor2vel_timer, vtend_timer
    use parcel_interpl, only : grid2par_timer, par2grid_timer
    use parcel_init, only : init_timer
    use ls_rk, only : rk_timer
    implicit none

    integer :: epic_timer

    contains

        subroutine register_all_timers
            call register_timer('epic', epic_timer)
            call register_timer('parcel container resize', resize_timer)
            call register_timer('par2grid', par2grid_timer)
            call register_timer('grid2par', grid2par_timer)
            call register_timer('parcel split', split_timer)
            call register_timer('parcel merge', merge_timer)
            call register_timer('laplace correction', lapl_corr_timer)
            call register_timer('gradient correction', grad_corr_timer)
            call register_timer('net vorticity correction', vort_corr_timer)
            call register_timer('parcel initialisation', init_timer)
            call register_timer('parcel diagnostics', parcel_stats_timer)
            call register_timer('parcel I/O', parcel_io_timer)
            call register_timer('parcel diagnostics I/O', parcel_stats_io_timer)
            call register_timer('field I/O', field_io_timer)
            call register_timer('field diagnostics', field_stats_timer)
            call register_timer('field diagnostics I/O', field_stats_io_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('parcel push', rk_timer)
            call register_timer('merge nearest', merge_nearest_timer)
            call register_timer('merge tree resolve', merge_tree_resolve_timer)
        end subroutine register_all_timers

end module test_utils
