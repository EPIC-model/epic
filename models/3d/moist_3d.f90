! =============================================================================
!                       MPIC plume test case
!
! This module enerates a spherical plume of approximately uniform liquid-wate
! buoyancy and specific humidity touching the lower surface, together with a
! stratified zone aloft.
! =============================================================================
module moist_3d
    use phys_constants
    use phys_parameters
    use constants
    use h5_writer
    implicit none

    private

    double precision, allocatable :: buoyg(:, :, :), humg(:, :, :)

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

    type(plume_type) :: moist

    public :: moist_init, &
              moist

    contains

        subroutine moist_init(h5handle, nx, ny, nz, origin, dx)
            integer(hid_t),   intent(inout) :: h5handle
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            integer                         :: i, j, k
            double precision                :: r2, pos(3), zbot
            double precision                :: rhb, bpl, dbdz, zb, hen, hpl, radsq

            if (moist%rhb .gt. one) then
                write(*,*) ' The relative humidity (as a fraction) must be < 1.'
                write(*,*)  '*** stopping ***'
                stop
            endif



            ! the specific humidity fraction in the plume:
            hpl = dexp(- moist%zc)

            write(*,'(a,f7.5)') &
                '  The specific humidity fraction inside the plume q_pl = ', hpl

            ! the specific humidity fraction outside the plume q_bg
            hen = moist%mu * hpl

            write(*,'(a,f7.5)') &
                '  The specific humidity fraction outside the plume q_bg = ', hen

            !From H and q_bg, obtain the height of the top of the mixed layer z_b:'
            zb = dlog(rhb / hen)

            write(*,'(a,f7.5)') &
                '  The height of the top of the mixed layer z_b = ln(H/q_bg) = ', zb

            ! From z_c, z_d & z_m, obtain the squared buoyancy frequency above z = z_b:
            dbdz = glat * (hpl - dexp(-moist%zm)) / (moist%zm - moist%zd)

            write(*,'(a,f7.5)') &
                '  The buoyancy frequency in the stratified zone is ', dsqrt(dbdz)

            !Also obtain the plume liquid-water buoyancy (using also z_b):
            bpl = dbdz * (moist%zd - zb)
            write(*,'(a,f7.5)') '  The plume liquid water buoyancy b_pl = ', bpl
            write(*,'(a,f7.5)') &
                '  corresponding to (theta_l-theta_l0)/theta_l0 = ', bpl / gravity


            if (two * moist%rad .gt. zb) then
                write(*,'(a,f7.3)') '  The radius must be less than z_b/2 = ', zb / two
                write(*,*)  '*** stopping ***'
                stop
            endif
            radsq = moist%rad ** 2

            ! Scale by radsq for use below:
            moist%e = moist%e / radsq

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(humg(0:nz, 0:ny-1, 0:nx-1))

            buoyg = zero
            humg = zero

            zbot = moist%zmin + zb

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz

                        pos = origin + dx * dble((/i, j, k/))

                        r2 = pos(1) ** 2 &
                           + pos(2) ** 2 &
                           + pos(3) ** 2

                        if (r2 <= radsq) then
                            buoyg(k, j, i) = bpl * (one + moist%e(1) * pos(1) * pos(2)  &
                                                        + moist%e(2) * pos(1) * pos(3)  &
                                                        + moist%e(3) * pos(2) * pos(3))
                            humg(k, j, i) = hpl
                        else
                            if (pos(3) < zb) then
                                ! Mixed layer:
                                buoyg(k, j, i) = zero
                                humg(k, j, i) = hen
                            else
                                ! Stratified layer
                                buoyg(k, j, i)= dbdz * (pos(3) - zbot)
                                humg(k, j, i) = rhb * dexp(moist%zmin - pos(3))
                            endif
                        endif
                    enddo
                enddo
            enddo

            call write_h5_dataset(h5handle, '/', 'buoyancy', buoyg)
            call write_h5_dataset(h5handle, '/', 'humidity', humg)

            deallocate(buoyg)
            deallocate(humg)

        end subroutine moist_init

end module moist_3d
