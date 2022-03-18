! =============================================================================
!                       MPIC plume test case
!
! This module enerates a spherical plume of approximately uniform liquid-wate
! buoyancy and specific humidity touching the lower surface, together with a
! stratified zone aloft.
! =============================================================================
module moist_3d
    use physics, only : write_physical_quantities, &
                        height_c, q_0, gravity, theta_0, L_v, c_p
    use constants
    use netcdf_writer
    implicit none

    private

    double precision, allocatable :: buoyg(:, :, :), humg(:, :, :)

    integer :: buo_id, hum_id

    type plume_type
        double precision :: H                   ! Relative humidity H (as a fraction < 1) at z_b
        double precision :: z_c                 ! Lifting condensation level
        double precision :: mu                  ! q_bg/q_pl where q_bg is the specific humidity fraction
                                                ! in the mixed layer outside of the plume (note, this should be
                                                ! between H and 1; a value < H would put z_c below z_b)
        double precision :: z_d                 ! Level of dry neutral stratification z_n
        double precision :: z_m                 ! Level of moist neutral stratification (cloud top) z_t
        double precision :: r_plume             ! Radius of the plume, R (must be < z_b/2)
        double precision :: e_values(3)         ! To create asymmetry, we vary the buoyancy in the plume
                                                ! according to  b = b_pl*[1 + (e1*x*y+e2*x*z+e3*yz)/R^2].
    end type plume_type

    type(plume_type) :: moist

    public :: moist_init, &
              moist

    contains

        subroutine moist_init(ncid, dimids, nx, ny, nz, origin, dx)
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            integer                         :: i, j, k
            double precision                :: rpos1, rpos2, rpos3, r2, pos(3), centre(3), extent(3)
            double precision                :: b_pl, dbdz, z_b, h_bg, h_pl, radsq

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='buoyancy',                &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=buo_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='humidity',                     &
                                       long_name='humidity',                &
                                       std_name='',                         &
                                       unit='-',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=hum_id)

            call close_definition(ncid)

            if (moist%H .gt. one) then
                write(*,*) 'Error: Relative humidity fraction must be < 1.'
                stop
            endif

            h_pl = q_0 * dexp(- moist%z_c / height_c)

            write(*,"('Humidity inside the plume is ',f6.3)") h_pl

            if ((moist%mu > one) .or. (moist%mu <= moist%H)) then
                write(*,'("Error: mu must be between H and 1. The selected value is ", f6.3)') moist%mu
                stop
            endif

            h_bg = moist%mu * h_pl

            write(*,"('Background humidity is ',f6.3)") h_bg

            z_b = height_c * dlog(q_0 * moist%H / h_bg)

            write(*,"('Base of mixed layer is ',f12.3)") z_b

            dbdz = (gravity * L_v / (c_p * theta_0)) &
                 * (h_pl - q_0 * dexp(-moist%z_m / height_c)) / (moist%z_m - moist%z_d)

            write(*,"('The buoyancy frequency in the stratified zone is ',f12.3)") dsqrt(dbdz)

            !Also obtain the plume liquid-water buoyancy (using also z_b):
            b_pl = dbdz * (moist%z_d - z_b)
            write(*,'(a,f7.5)') '  The plume liquid water buoyancy b_pl = ', b_pl
            write(*,'(a,f7.5)') '  corresponding to (theta_l-theta_0)/theta_0 = ', b_pl * gravity


            if (two * moist%r_plume > z_b) then
                write(*,"('Error: Plume radius is too big. At most it can be ',f7.5)") f12 * z_b
                stop
            endif

            radsq = moist%r_plume ** 2

            moist%e_values = moist%e_values / radsq

            if (origin(3) > zero) then
                write(*, *) "The vertical origin must be zero."
                stop
            endif

            extent = dx * dble((/nx, ny, nz/))

            write(*,*) "Box layout:"
            write(*,*) "z_max           =", extent(3)
            write(*,*) "z_m             =", moist%z_m
            write(*,*) "z_d             =", moist%z_d
            write(*,*) "z_c             =", moist%z_c
            write(*,*) "z_b             =", z_b
            write(*,*) "zmin            =", origin(3)
            write(*,*) "top of plume    =", two * moist%r_plume
            write(*,*) "bottom of plume =", zero

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(humg(0:nz, 0:ny-1, 0:nx-1))

            buoyg = zero
            humg = zero

            centre = f12 * (two * origin + extent)

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz

                        pos = origin + dx * dble((/i, j, k/))

                        rpos1 = (pos(1) - centre(1))
                        rpos2 = (pos(2) - centre(2))
                        rpos3 = (pos(3) - moist%r_plume)
                        r2 = rpos1 ** 2 + rpos2 ** 2 + rpos3 ** 2

                        if (r2 <= radsq) then
                            buoyg(k, j, i) = b_pl * (one + moist%e_values(1) * rpos1 * rpos2  &
                                                         + moist%e_values(2) * rpos1 * rpos3  &
                                                         + moist%e_values(3) * rpos2 * rpos3)
                            humg(k, j, i) = h_pl
                        else
                            if (pos(3) < z_b) then
                                ! Mixed layer:
                                buoyg(k, j, i) = zero
                                humg(k, j, i) = h_bg
                            else
                                ! Stratified layer
                                buoyg(k, j, i)= dbdz * (pos(3) - z_b)
                                humg(k, j, i) = q_0 * moist%H * dexp(- pos(3) / height_c)
                            endif
                        endif
                    enddo
                enddo
            enddo

            call write_netcdf_dataset(ncid, buo_id, buoyg)
            call write_netcdf_dataset(ncid, hum_id, humg)

            call write_physical_quantities(ncid)

            deallocate(buoyg)
            deallocate(humg)

        end subroutine moist_init

end module moist_3d
