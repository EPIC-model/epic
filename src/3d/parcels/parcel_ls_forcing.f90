! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_ls_forcing
    use constants, only : zero, one, pi
    use mpi_timer, only : start_timer, stop_timer
    use options, only : parcel
    use parameters, only : nz, dxi
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid
    use fields
    use omp_lib
    use mpi_utils, only : mpi_exit_on_error
    use parcel_interpl, only : saturation_adjustment, par2grid_diag, trilinear

    implicit none

    private

    contains

        subroutine ls_tendencies(dt)
            double precision, intent(in)  :: dt
           call apply_subsidence_and_vorticity_adjustment(dt)
           call apply_ls_tendencies(dt)
           call saturation adjustment
        end subroutine ls_tendencies
        !
        ! @pre
        subroutine apply_subsidence_and_vorticity_adjustment(dt)
            double precision, intent(in)  :: dt
            double precision              :: weights(0:1,0:1,0:1)
            double precision              :: xf, yf, zf, xfc, yfc, zfc
            double precision              :: thetagradz, qcgradz, qlgradz, ww
            integer                       :: n, nc, is, js, ks
            double precision              :: xvortbar(0:nz), yvortbar(0:nz)

            ! use par2grid_diag theta
            ! use par2grid_diag ql+qv
            ! to be replaced by special par2grid?
            call par2grid_diag(.false.)

            ! code to apply subsidence, use local gradients

            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights, xf, yf, zf, xfc, yfc, zfc, thetagradz, qvgradz, qlgradz, ww)
            do n = 1, n_parcels

                call trilinear(parcels%position(:, n), is, js, ks, weights)

                xf = weights(0,0,1) + weights(0,1,1) + weights(1,0,1) + weights(1,1,1) ! fractional position along x
                yf = weights(0,1,0) + weights(0,1,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along y
                zf = weights(1,0,0) + weights(1,0,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along z
                xfc=one-xf
                yfc=one-yf
                zfc=one-zf

                thetagradz = dxi(3)*(&
                            xfc * (&
                                    yfc * (thetag(ks+1, js, is  ) - thetag(ks, js,   is))  &
                             +      yf  * (thetag(ks+1, js+1, is) - thetag(ks, js+1, is))) &
                     +      xf *  (&
                                    yfc * (thetag(ks+1, js,   is+1) - thetag(ks, js,   is+1))  &
                             +      yf  * (thetag(ks+1, js+1, is+1) - thetag(ks, js+1, is+1))))

                qvgradz = dxi(3)*(&
                            xfc * (&
                                    yfc * (qvg(ks+1, js, is  ) - qvg(ks, js,   is))  &
                             +      yf  * (qvg(ks+1, js+1, is) - qvg(ks, js+1, is))) &
                     +      xf *  (&
                                    yfc * (qvg(ks+1, js,   is+1) - qvg(ks, js,   is+1))  &
                             +      yf  * (qvg(ks+1, js+1, is+1) - qvg(ks, js+1, is+1))))

                qlgradz = dxi(3)*(&
                            xfc * (&
                                    yfc * (qlg(ks+1, js, is  ) - qlg(ks, js,   is))  &
                             +      yf  * (qlg(ks+1, js+1, is) - qlg(ks, js+1, is))) &
                     +      xf *  (&
                                    yfc * (qlg(ks+1, js,   is+1) - qlg(ks, js,   is+1))  &
                             +      yf  * (qlg(ks+1, js+1, is+1) - qlg(ks, js+1, is+1))))

                ww=get_bomex_subsidence_velocity(parcels%position(3, n))

                parcels%theta(n) = parcels%theta(n)-ww*thetagradz
                parcels%qv(n) = parcels%qv(n)-ww*qvgradz-ww*qlgradz
            enddo
            !$omp end do
            !$omp end parallel

            ! code to apply vorticity adjustment
            do nc = 1, 2
                call field_decompose_physical(vortg(:, :, :, nc), svor(:, :, :, nc))
            enddo

            ! ubar and vbar are used here to store the mean x and y components of the vorticity
            ! may need a correction factor
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                xvortbar = svor(:, 0, 0, 1)
                yvortbar = svor(:, 0, 0, 2)
            endif

            !$omp parallel default(shared)parcels%position(:, n)
            !$omp do private(n, is, js, ks, weights, zf, zfc)
            do n = 1, n_parcels

                call trilinear(parcels%position(:, n), is, js, ks, weights)

                zf = weights(1,0,0) + weights(1,0,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along z
                zfc= one-zf

                parcels%vorticity(1, n) = parcels%vorticity(1, n) - zfc*xvortbar(ks) - zf*xvortbar(ks+1) + get_bomex_xvort(parcels%position(3, n))
                parcels%vorticity(2, n) = parcels%vorticity(2, n) - zfc*yvortbar(ks) - zf*yvortbar(ks+1) + get_bomex_yvort(parcels%position(3, n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_subsidence_and_vorticity_adjustment

        !
        ! @pre
        subroutine apply_ls_tendencies
            double precision, intent(in)  :: dt

            !$omp parallel default(shared)parcels%position(:, n)
            !$omp do private(n,)
            do n = 1, n_parcels
                parcels%theta(n) = parcels%theta(n) + get_bomex_theta_ls(parcels%position(3, n))*dt
                parcels%qv(n) = parcels%qv(n) + get_bomex_qt_ls(parcels%position(3, n))*dt
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_ls_tendencies
        !
        ! @pre
        pure function get_bomex_subsidence_velocity(zz):
           if (zz < 1500.0) then
              ww = zz*(-0.0065)/1500.0
           else if (z[k] < 2100.0) then
              ww = -0.0065 + (zz-1500.0)*(0.0065)/(2100.0-1500.0)
           else
              ww = 0.0
           endif
           return ww
        end function get_bomex_subsidence_velocity

        !
        ! @pre
        pure function get_bomex_theta_ls(zz)
          double precision, intent(in) :: zz
          double precision, intent(out) :: theta_ls
          double precision, parameter :: iday=1./86400.

          if (zz < 1500.0) then
             theta_ls = -2.0*iday
          else
             theta_ls = (-2.0 + (zz-1500.0)*(2.0/(3000.0-1500.0)))*iday
          endif
          return theta_ls
        end function get_bomex_theta_ls

        !
        ! @pre
        elemental function get_bomex_qt_ls(zz) result qt_ls
          double precision, intent(in) :: zz
          double precision, intent(out) :: qt_ls

          if(zz <= 300.0) then
            qt_ls = -1.2e-8
          else if(zz <= 500.0) then
            qt_ls = -1.2e-8 + (zz-300.)*(1.2e-8)/(500.0-300.0)
          else
            qt_ls = 0.0
          endif
          return qt_ls
        end function get_bomex_qt_ls

        !
        ! @pre
        elemental function  get_bomex_xvort(zz) result xvort
          double precision, intent(in) :: zz
          double precision, intent(out) :: xvort

          if (zz<300.0) then
              xvort=0.0043*cos(pi*(zz+600)/1200.0)+0.0023
          else if (zz<600.0) then
              xvort=0.0043*cos(pi*(zz+600)/1200.0)+0.0023
          else if (zz<1400.0) then
              xvort=-0.001-0.001*cos(pi*(zz-600)/800.0)
          else
              xvort=0.0
          end if

        end function get_bomex_xvort

        ! @pre
        elemental function get_bomex_yvort(zz) result yvort
          double precision, intent(in) :: zz
          double precision, intent(out) :: yvort

          if(zz<300.0) then
              yvort=0.0062*cos((pi*(zz-300.0)/600.0)**2)-0.0067
          else if (zz<600.0) then
              yvort=0.0027*cos(pi*(zz-300.0)/600.0)-0.0032
          else if (zz<1400.0) then
              yvort=0.005*cos(pi*(zz-1400.0)/1600.0)-0.0032
          else
              yvort=0.0018
          end if

        end subroutine get_bomex_yvort

end module parcel_ls_forcing
