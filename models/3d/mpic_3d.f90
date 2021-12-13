! =============================================================================
!                       MPIC plume test case
!
! This module enerates a spherical plume of approximately uniform liquid-wate
! buoyancy and specific humidity touching the lower surface, together with a
! stratified zone aloft.
! =============================================================================



! read(*,*) mu
! hen=mu*hpl
! write(*,'(a,f7.5)') &
!   '  The specific humidity fraction outside the plume q_bg = ',hen
!
!  !From H and q_bg, obtain the height of the top of the mixed layer z_b:'
! zb=log(rhb/hen)
! write(*,'(a,f7.5)') &
!   '  The height of the top of the mixed layer z_b = ln(H/q_bg) = ',zb
!
! write(*,*)
! write(*,*) ' Level of dry neutral stratification z_n?'
! read(*,*) zd
! write(*,*) ' Level of moist neutral stratification (cloud top) z_t?'
! read(*,*) zm
!
!  !From z_c, z_d & z_m, obtain the squared buoyancy frequency above z = z_b:
! dbdz=glat*(hpl-exp(-zm))/(zm-zd)
! write(*,'(a,f7.5)') '  The buoyancy frequency in the stratified zone is ', &
!      sqrt(dbdz)
!
!  !Also obtain the plume liquid-water buoyancy (using also z_b):
! bpl=dbdz*(zd-zb)
! write(*,'(a,f7.5)') '  The plume liquid water buoyancy b_pl = ',bpl
! write(*,'(a,f7.5)') '  corresponding to (theta_l-theta_l0)/theta_l0 = ', &
!      bpl/gravity
!
!  !Specify plume size:
! write(*,*)
! write(*,*) ' Radius of the plume, R (must be < z_b/2)?'
! read(*,*) rad
! if (two*rad .gt. zb) then
!   write(*,'(a,f7.3)') '  The radius must be less than z_b/2 = ',zb/two
!   write(*,*)  '*** stopping ***'
!   stop
! endif
! radsq=rad**2
!
! write(*,*)
! write(*,*) ' To create asymmetry, we vary the buoyancy in the plume'
! write(*,*) ' according to  b = bpl*[1 + (e1*x*y+e2*x*z+e3*yz)/R^2].'
! write(*,*) ' Enter e1, e2 and e3:'
! read(*,*) e1,e2,e3
!  !Scale by radsq for use below:
! e1=e1/radsq
! e2=e2/radsq
! e3=e3/radsq
!
! write(*,*) ' The initial flow is (u,v,w) = (S*z,0,0); enter the shear S:'
! read(*,*) shear


module mpic_3d
    use phys_constants
    use constants
    use h5_writer
    implicit none

    private

    double precision, allocatable :: buoyg(:, :, :), humg(:, :, :)

    double precision :: ref_theta = 303.15d0    ![K] reference potential temperature

    type plume_type
        double precision :: rhb   ! Relative humidity H (as a fraction < 1) at z_b
        double precision :: zc    ! Lifting condensation level
        double precision :: mu    ! q_bg/q_pl where q_bg is the specific humidity fraction
                                  ! in the mixed layer outside of the plume (note, this should be
                                  ! between H and 1; a value < H would put z_c below z_b)
        double precision :: zd    ! Level of dry neutral stratification z_n
        double precision :: zm    ! Level of moist neutral stratification (cloud top) z_t
        double precision :: rad   ! Radius of the plume, R (must be < z_b/2)
        double precision :: e(3)  ! To create asymmetry, we vary the buoyancy in the plume
                                  ! according to  b = bpl*[1 + (e1*x*y+e2*x*z+e3*yz)/R^2].
        double precision :: shear ! The initial flow is (u,v,w) = (S*z,0,0); enter the shear S
        double precision :: zmin
    end type plume_type

    type(flow_type) :: mpic_flow

    public :: mpic_init, &
              mpic_flow

    contains

        subroutine mpic_init(h5handle, nx, ny, nz, origin, dx)
            integer(hid_t),   intent(inout) :: h5handle
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            integer                         :: k
            type(bubble_type)               :: bubble

            if (mpic_flow%rhb .gt. one) then
                write(*,*) ' The relative humidity (as a fraction) must be < 1.'
                write(*,*)  '*** stopping ***'
                stop
            endif



            ! the specific humidity fraction in the plume:
            hpl = dexp(- mpic_flow%zc)

            write(*,'(a,f7.5)') &
                '  The specific humidity fraction inside the plume q_pl = ', hpl

            ! the specific humidity fraction outside the plume q_bg
            hen = mpic_flow%mu * hpl

            write(*,'(a,f7.5)') &
                '  The specific humidity fraction outside the plume q_bg = ', hen

            !From H and q_bg, obtain the height of the top of the mixed layer z_b:'
            zb = dlog(rhb / hen)

            write(*,'(a,f7.5)') &
                '  The height of the top of the mixed layer z_b = ln(H/q_bg) = ', zb

            ! From z_c, z_d & z_m, obtain the squared buoyancy frequency above z = z_b:
            dbdz = glat * (hpl - dexp(-zm)) / (zm - zd)

            write(*,'(a,f7.5)') &
                '  The buoyancy frequency in the stratified zone is ', dsqrt(dbdz)

            !Also obtain the plume liquid-water buoyancy (using also z_b):
            bpl = dbdz * (zd - zb)
            write(*,'(a,f7.5)') '  The plume liquid water buoyancy b_pl = ', bpl
            write(*,'(a,f7.5)') &
                '  corresponding to (theta_l-theta_l0)/theta_l0 = ', bpl / gravity


            if (two * rad .gt. zb) then
                write(*,'(a,f7.3)') '  The radius must be less than z_b/2 = ', zb / two
                write(*,*)  '*** stopping ***'
                stop
            endif
            radsq = rad ** 2

            ! Scale by radsq for use below:
            e = e / radsq




            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(humg(0:nz, 0:ny-1, 0:nx-1))

            buoyg = zero
            humg = zero

            call mpic_init(nx, ny, nz, origin, dx)

            call write_h5_dataset(h5handle, '/', 'buoyancy', buoyg)
            call write_h5_dataset(h5handle, '/', 'humidity', humg)

            deallocate(buoyg)
            deallocate(humg)

        end subroutine mpic_init

        subroutine mpic_init(nx, ny, nz, origin, dx)
            integer,           intent(in) :: nx, ny, nz
            double precision,  intent(in) :: origin(3), dx(3)
            integer                       :: i, j, k
            double precision              :: r2, pos(3), zbot

            zbot = zmin + zb

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz

                        pos = origin + dx * dble((/i, j, k/))

                        r2 = pos(1) ** 2 &
                           + pos(2) ** 2 &
                           + pos(3) ** 2

                        if (r2 <= radsq) then
                            buoyg(k, j, i) = bpl * (one + e(1) * pos(1) * pos(2)  &
                                                        + e(2) * pos(1) * pos(3)  &
                                                        + e(3) * pos(2) * pos(3))
                            hmg(k, j, i) = hpl
                        else
                            if (pos(3) < zb) then
                                ! Mixed layer:
                                buoyg(k, j, i) = zero
                                humg(k, j, i) = hen
                            else
                                ! Stratified layer
                                buoyg(k, j, i)= dbdz * (pos(3) - zbot)
                                humg(k, j, i) = rhb * dexp(zmin - pos(3))
                            endif
                        endif
                    enddo
                enddo
            enddo
        end subroutine mpic_init
end module mpic_3d
