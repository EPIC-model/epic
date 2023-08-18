! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_ls_forcings
    use constants, only : zero, one, pi
    use mpi_timer, only : start_timer, stop_timer
    use options, only : parcel
    use parameters, only : nx,ny,nz, dxi, dx
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid
    use fields
    use omp_lib
    use inversion_utils, only : field_decompose_physical
    use mpi_utils, only : mpi_exit_on_error
    use mpi_layout
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce 
    use parcel_interpl, only : saturation_adjustment, par2grid, par2grid_diag, trilinear

    implicit none

    private

    public :: apply_ls_forcings

    contains

        subroutine apply_ls_forcings(dt)
            double precision, intent(in)  :: dt
           call apply_subsidence_and_vorticity_adjustment(dt)
           call apply_ls_tendencies(dt)
           call saturation_adjustment
        end subroutine apply_ls_forcings
        !
        ! @pre
        subroutine apply_subsidence_and_vorticity_adjustment(dt)
            double precision, intent(in)  :: dt
            double precision              :: weights(0:1,0:1,0:1)
            double precision              :: xf, yf, zf, xfc, yfc, zfc
            double precision              :: thetagradz, qvgradz, qlgradz, ww
            integer                       :: n, kk, is, js, ks
            double precision              :: xvortbar(0:nz), yvortbar(0:nz), nudgefac

            ! use par2grid_diag theta
            ! use par2grid_diag ql+qv
            ! to be replaced by special par2grid?
            call par2grid(.false.)
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

                parcels%theta(n) = parcels%theta(n)-dt*ww*thetagradz
                parcels%qv(n) = parcels%qv(n)-dt*ww*(qvgradz+qlgradz)
            enddo
            !$omp end do
            !$omp end parallel

            do kk=0,nz
                xvortbar(kk)= sum(vortg(kk,box%lo(2):box%hi(2),box%lo(1):box%hi(1),1))
                yvortbar(kk)= sum(vortg(kk,box%lo(2):box%hi(2),box%lo(1):box%hi(1),2))
            end do

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               xvortbar,                &
                               nz,                      &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               yvortbar,                &
                               nz,                      &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            do kk=0,nz
               xvortbar(kk)= xvortbar(kk)/(nx*ny)-get_bomex_xvort(lower(3)+kk*dx(3))
               yvortbar(kk)= yvortbar(kk)/(nx*ny)-get_bomex_yvort(lower(3)+kk*dx(3))
            enddo

            !if(cart%rank == cart%root) then
            !    write(*,*) 'xvortbar'
            !    write(*,*) xvortbar
            !    write(*,*) 'yvortbar'
            !    write(*,*) yvortbar
            !end if

            ! nudge velocity profiles with 120 seconds time scale 
            nudgefac=(1.0-exp(-dt/120.0))

            !$omp parallel default(shared)
            !$omp do private(n, is, js, ks, weights, zf, zfc)
            do n = 1, n_parcels
            
                call trilinear(parcels%position(:, n), is, js, ks, weights)

                zf = weights(1,0,0) + weights(1,0,1) + weights(1,1,0) + weights(1,1,1) ! fractional position along z
                zfc= one-zf

                parcels%vorticity(1, n) = parcels%vorticity(1, n) - zfc*xvortbar(ks)*nudgefac &
                                          - zf*xvortbar(ks+1)*nudgefac
                parcels%vorticity(2, n) = parcels%vorticity(2, n) - zfc*yvortbar(ks)*nudgefac &
                                          - zf*yvortbar(ks+1)*nudgefac
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_subsidence_and_vorticity_adjustment

        !
        ! @pre
        subroutine apply_ls_tendencies(dt)
            double precision, intent(in)  :: dt
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%theta(n) = parcels%theta(n) + get_bomex_theta_ls(parcels%position(3, n))*dt
                parcels%qv(n) = parcels%qv(n) + get_bomex_qt_ls(parcels%position(3, n))*dt
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_ls_tendencies
        !
        ! @pre
        elemental function get_bomex_subsidence_velocity(zz) result (ww)
          double precision, intent(in) :: zz
          double precision :: ww

          if (zz < 1500.0) then
             ww = zz*(-0.0065)/1500.0
          else if (zz < 2100.0) then
             ww = -0.0065 + (zz-1500.0)*(0.0065)/(2100.0-1500.0)
          else
             ww = 0.0
          endif
        end function get_bomex_subsidence_velocity

        !
        ! @pre
        elemental function get_bomex_theta_ls(zz) result (theta_ls)
          double precision, intent(in) :: zz
          double precision :: theta_ls
          double precision, parameter :: iday=1./86400.

          if (zz < 1500.0) then
             theta_ls = -2.0*iday
          else
             theta_ls = (-2.0 + (zz-1500.0)*(2.0/(3000.0-1500.0)))*iday
          endif
        end function get_bomex_theta_ls

        !
        ! @pre
        elemental function get_bomex_qt_ls(zz) result (qt_ls)
          double precision, intent(in) :: zz
          double precision :: qt_ls

          if(zz <= 300.0) then
            qt_ls = -1.2e-8
          else if(zz <= 500.0) then
            qt_ls = -1.2e-8 + (zz-300.)*(1.2e-8)/(500.0-300.0)
          else
            qt_ls = 0.0
          endif
        end function get_bomex_qt_ls

        !
        ! @pre
        elemental function  get_bomex_xvort(zz) result (xvort)
          double precision, intent(in) :: zz
          double precision :: xvort

          if (zz<600.0) then
              xvort=0.0043*cos(pi*(zz+600)/1200.0)+0.0023
          else if (zz<1400.0) then
              xvort=-0.001-0.001*cos(pi*(zz-600)/800.0)
          else
              xvort=0.0
          end if

        end function get_bomex_xvort

        ! @pre
        elemental function get_bomex_yvort(zz) result (yvort)
          double precision, intent(in) :: zz
          double precision :: yvort

          if(zz<300.0) then
              yvort=0.0062*cos((pi*(zz-300.0)/600.0)**2)-0.0067
          else if (zz<600.0) then
              yvort=0.0027*cos(pi*(zz-300.0)/600.0)-0.0032
          else if (zz<1400.0) then
              yvort=0.005*cos(pi*(zz-1400.0)/1600.0)-0.0032
          else
              yvort=0.0018
          end if

        end function get_bomex_yvort

end module parcel_ls_forcings
