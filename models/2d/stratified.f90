! =============================================================================
!                       Just set up a linear stratification
!
! =============================================================================

module stratified
    use constants
    use netcdf_writer
    use physics, only : write_physical_quantities, &
                        p_surf, pressure_scale_height, p_ref, r_d, c_p

    implicit none
    double precision, parameter :: tk0c = 273.15       ! Temperature of freezing in Kelvin
    double precision, parameter :: qsa1 = 3.8          ! Top in equation to calculate qsat
    double precision, parameter :: qsa2 = -17.2693882  ! Constant in qsat equation
    double precision, parameter :: qsa3 = 35.86        ! Constant in qsat equation
    double precision, parameter :: qsa4 = 6.109        ! Constant in qsat equation

    private

    double precision, allocatable :: thetag(:, :), qvg(:, :)

    integer :: theta_id, qv_id

    public :: stratified_init

    contains

        subroutine stratified_init(ncid, dimids, nx, nz, origin, dx)
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, nz
            double precision, intent(in)    :: origin(2)
            double precision, intent(in)    :: dx(2)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='theta',                     &
                                       long_name='theta',                &
                                       std_name='',                         &
                                       unit='K',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=theta_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='qv',                     &
                                       long_name='water vapour mixing ratio',                &
                                       std_name='',                         &
                                       unit='kg/kg',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=qv_id)

            call close_definition(ncid)

            allocate(thetag(0:nz, 0:nx-1))
            allocate(qvg(0:nz, 0:nx-1))

            call stratified_field_init(nx, nz, origin, dx)

            call write_netcdf_dataset(ncid, theta_id, thetag)
            call write_netcdf_dataset(ncid, qv_id, qvg)

            call write_physical_quantities(ncid)

            deallocate(thetag)
            deallocate(qvg)

        end subroutine stratified_init

       subroutine stratified_field_init(nx, nz, origin, dx)
            integer,           intent(in) :: nx, nz
            double precision,  intent(in) :: origin(2), dx(2)
            double precision              :: pos(2)
            integer                       :: i, j
            double precision :: rh, temp, exn, press, qsat, theta, qv, theta_surf

            exn=(p_surf/p_ref)**(r_d/c_p)
            theta_surf=283.0/exn
            rh=0.50

            do j = 0, nz
                do i = 0, nx-1
                    pos = origin + dx * dble((/i, j/))
                    theta = theta_surf+1e-3*pos(2)
                    press=p_surf*exp(-pos(2)/pressure_scale_height)
                    exn=(press/p_ref)**(r_d/c_p)
                    temp=theta*exn
                    qsat = qsa1/(0.01*press*exp(qsa2*(temp - tk0c)/(temp - qsa3)) - qsa4)
                    qv = rh*qsat
                    thetag(j ,i) = theta
                    qvg(j, i) = qv
                enddo
            enddo

        end subroutine stratified_field_init

end module
