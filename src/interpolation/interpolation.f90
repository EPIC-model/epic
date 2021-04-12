! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module interpolation
    use constants, only : max_num_parcels
    use parameters, only : nx, nz
    use options, only : parcel_info, interpl
    use parcel_container, only : parcel_container_type, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use ellipse
    use fields
    use taylorgreen, only : get_flow_velocity, &
                            get_flow_gradient

    implicit none

    private :: par2grid_elliptic,                       &
               par2grid_non_elliptic,                   &
               grid2par_elliptic,                       &
               grid2par_non_elliptic,                   &
               cache_parcel_interp_weights_elliptic,    &
               cache_parcel_interp_weights_non_elliptic


    ! flag to see if this is intiialised or not
    logical :: initialised = .false.

    ! cached variables
    integer, allocatable, dimension(:, :) ::  is, js

    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer ij(2, 4)

    ! trilinear interpolation requires 4 grid points in 2D
    integer, parameter :: ngp = 4

    ! interpolation weights
    double precision, allocatable, dimension(:, :) :: weights

    private :: is, js, weights, ngp

    contains

        subroutine initialise_parcel_interp
            if (initialised) then
                return
            endif

            initialised = .true.

            ! allocate cache arrays
            if (parcel_info%is_elliptic) then
                ! we have 2 points per parcel
                allocate(is(ngp, 2 * max_num_parcels))
                allocate(js(ngp, 2 * max_num_parcels))
                allocate(weights(ngp, 2 * max_num_parcels))
            else
                allocate(is(ngp, max_num_parcels))
                allocate(js(ngp, max_num_parcels))
                allocate(weights(ngp, max_num_parcels))
            endif

        end subroutine initialise_parcel_interp

        ! finalise parcel interpolation data structures
        subroutine finalise_parcel_interp
            if (.not. initialised) then
                return
            endif

            deallocate(is)
            deallocate(js)
            deallocate(weights)

        end subroutine finalise_parcel_interp

        subroutine par2grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(0:, -1:, :)

            field = 0.0

            if (parcel_info%is_elliptic) then
                call par2grid_elliptic(parcels, attrib, field)
            else
                call par2grid_non_elliptic(parcels, attrib, field)
            endif

            ! apply free slip boundary condition
            field(:, 0, :)  = 2.0 * field(:, 0, :)
            field(:, nz, :) = 2.0 * field(:, nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            field(:, 1, :)    = field(:, 1, :) + field(:, -1, :)
            field(:, nz-1, :) = field(:, nz-1, :) + field(:, nz+1, :)

        end subroutine par2grid

        subroutine par2grid_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(0:, -1:, :)
            integer                                    :: ncomp, ngp
            double precision                           :: points(2, 2)
            integer                                    :: n, p, c, i, k
            integer                                    :: the_shape(3)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n, 1),   &
                                            parcels%B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2

                    k = 2 * (n - 1) + p

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            field(is(i, k), js(i, k), c) = field(is(i, k), js(i, k), c)     &
                                                         + weights(i, k) * attrib(n, c)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine par2grid_elliptic


        subroutine par2grid_non_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(0:, -1:, :)
            integer                                    :: ncomp, ngp
            integer                                    :: n, c, i
            integer                                    :: the_shape(3)
            double precision                           :: pos(2)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels
                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        field(is(i, n), js(i, n), c) = field(is(i, n), js(i, n), c) &
                                                     + weights(i, n) * attrib(n, c)
                    enddo
                enddo
            enddo

        end subroutine par2grid_non_elliptic


        subroutine grid2par(position, volume, attrib, field, B, exact)
            double precision,           intent(in)  :: position(:, :)
            double precision,           intent(in)  :: volume(:, :)
            double precision,           intent(out) :: attrib(:, :)
            double precision,           intent(in)  :: field(0:, -1:, :)
            double precision, optional, intent(in)  :: B(:, :)
            character(*), optional, intent(in)      :: exact

            if (parcel_info%is_elliptic) then
                if (.not. present(B)) then
                    print *, "B matrix not passed to grid2par!"
                    stop
                endif
                if(present(exact)) then
                   call grid2par_elliptic(position, volume, B, attrib, field, exact=exact)
                else
                   call grid2par_elliptic(position, volume, B, attrib, field)
                endif
            else
                if(present(exact)) then
                   call grid2par_non_elliptic(position, attrib, field, exact=exact)
                else
                   call grid2par_non_elliptic(position, attrib, field)
                endif
            endif

        end subroutine grid2par


        subroutine grid2par_add(position, volume, attrib, field, B, exact)
            double precision,           intent(in)  :: position(:, :)
            double precision,           intent(in)  :: volume(:, :)
            double precision,           intent(out) :: attrib(:, :)
            double precision,           intent(in)  :: field(0:, 0:, :)
            double precision, optional, intent(in)  :: B(:, :)
            character(*), optional, intent(in)      :: exact

            if (parcel_info%is_elliptic) then
                if (.not. present(B)) then
                    print *, "B matrix not passed to grid2par!"
                    stop
                endif
                if(present(exact)) then
                   call grid2par_elliptic(position, volume, B, attrib, field, add=.true., exact=exact)
                else
                   call grid2par_elliptic(position, volume, B, attrib, field, add=.true.)
                endif
            else
                if(present(exact)) then
                   call grid2par_non_elliptic(position, attrib, field, add=.true., exact=exact)
                else
                   call grid2par_non_elliptic(position, attrib, field, add=.true.)
                endif
            endif

        end subroutine grid2par_add


        subroutine grid2par_elliptic(position, volume, B, attrib, field, add, exact)
            double precision, intent(in)  :: position(:, :)
            double precision, intent(in)  :: volume(:, :)
            double precision, intent(in)  :: B(:, :)
            double precision, intent(out) :: attrib(:, :)
            double precision, intent(in)  :: field(0:, -1:, :)
            logical, optional, intent(in) :: add
            character(*), optional, intent(in)      :: exact
            integer                       :: ncomp, ngp
            double precision              :: points(2, 2)
            integer                       :: n, p, c, i, k
            integer                       :: the_shape(3)


            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                   attrib(1:n_parcels, :) = 0.0
               endif
            else
               attrib(1:n_parcels, :) = 0.0
            endif

            ! put if statement here for computational efficiency
            if(present(exact)) then
               do n = 1, n_parcels

                  points = get_ellipse_points(position(n, :), volume(n, 1), B(n, :))

                  do p = 1, 2
                     call apply_periodic_bc(points(p, :))
                     if(exact=='velocity') then
                        attrib(n,:)=attrib(n,:)+0.5*get_flow_velocity(points(p, :))
                     elseif(exact=='strain') then
                        attrib(n,:)=attrib(n,:)+0.5*get_flow_gradient(points(p, :))
                     else
                        print *, "Exact interpolation field passed not implemented"
                        stop
                     end if
                  end do
               end do
               return
            endif

            do n = 1, n_parcels
                ! we have 2 points per ellipse
                do p = 1, 2

                    k = 2 * (n - 1) + p

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            attrib(n, c) = attrib(n, c) &
                                         + weights(i, k) * field(is(i, k), js(i, k), c)
                        enddo
                    enddo
                enddo
            enddo

        end subroutine grid2par_elliptic


        subroutine grid2par_non_elliptic(position, attrib, field, add, exact)
            double precision, intent(in)  :: position(:, :)
            double precision, intent(out) :: attrib(:, :)
            double precision, intent(in)  :: field(0:, -1:, :)
            logical, optional, intent(in) :: add
            character(*), optional, intent(in)      :: exact
            integer                       :: ncomp, ngp
            integer                       :: n, c, i
            integer                       :: the_shape(3)
            double precision              :: pos(2)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                   attrib(1:n_parcels, :) = 0.0
               endif
            else
               attrib(1:n_parcels, :) = 0.0
            endif

            ! put if statement here for computational efficiency
            if(present(exact)) then
               do n = 1, n_parcels
                  pos = position(n, :)
                  call apply_periodic_bc(pos)
                  if(exact=='velocity') then
                     attrib(n,:)=attrib(n,:)+get_flow_velocity(pos)
                  elseif(exact=='strain') then
                     attrib(n,:)=attrib(n,:)+get_flow_gradient(pos)
                  else
                     print *, "Exact interpolation field passed not implemented"
                     stop
                  end if
               end do
               return
            endif

            do n = 1, n_parcels
                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        attrib(n, c) = attrib(n, c) &
                                     + weights(i, n) * field(is(i, n), js(i, n), c)
                    enddo
                enddo
            enddo

        end subroutine grid2par_non_elliptic

        ! cache indices and weights of tri-linear interpolation
        subroutine cache_parcel_interp_weights(parcels)
            type(parcel_container_type), intent(in) :: parcels

            if (parcel_info%is_elliptic) then
                call cache_parcel_interp_weights_elliptic(parcels)
            else
                call cache_parcel_interp_weights_non_elliptic(parcels)
            endif
        end subroutine cache_parcel_interp_weights

        ! cache indices and weights of tri-linear interpolation (elliptic version)
        subroutine cache_parcel_interp_weights_elliptic(parcels)
            type(parcel_container_type), intent(in) :: parcels
            double precision                        :: pos(2), xy(2)
            integer                                 :: idx(2), n, p
            double precision                        :: points(2, 2)

            do n = 1, n_parcels

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n, 1),   &
                                            parcels%B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2
                    ! the weight is halved due to 2 points per ellipse --> factor = 0.5
                    call trilinear_interp(points(p, :), 2 * (n - 1) + p, 0.5d0)
                enddo
            enddo
        end subroutine cache_parcel_interp_weights_elliptic

        ! cache indices and weights of tri-linear interpolation (non elliptic version)
        subroutine cache_parcel_interp_weights_non_elliptic(parcels)
            type(parcel_container_type), intent(in) :: parcels
            double precision                        :: pos(2), xy(2)
            integer                                 :: idx(2), n

            do n = 1, n_parcels
                call trilinear_interp(parcels%position(n, :), n, 1.0d0)
            enddo
        end subroutine cache_parcel_interp_weights_non_elliptic


        subroutine trilinear_interp(inpos, n, factor)
            double precision, intent(in) :: inpos(2), factor
            integer,          intent(in) :: n
            double precision             :: xy(2), pos(2)
            integer                      :: idx(2)

            pos = inpos

            ! ensure point is within the domain
            call apply_periodic_bc(pos)

            ! (i, j)
            idx = get_index(pos)
            xy = get_position(idx)
            weights(1, n) = factor * product(1.0 - abs(pos - xy) * dxi)
            is(1, n) = idx(1)
            js(1, n) = idx(2)

            ! (i+1, j)
            xy = get_position(idx + (/1, 0/))
            weights(2, n) = factor * product(1.0 - abs(pos - xy) * dxi)
            is(2, n) = idx(1) + 1
            js(2, n) = idx(2)

            ! (i, j+1)
            xy = get_position(idx + (/0, 1/))
            weights(3, n) = factor * product(1.0 - abs(pos - xy) * dxi)
            is(3, n) = idx(1)
            js(3, n) = idx(2) + 1

            ! (i+1, j+1)
            xy = get_position(idx + (/1, 1/))
            weights(4, n) = factor * product(1.0 - abs(pos - xy) * dxi)
            is(4, n) = idx(1) + 1
            js(4, n) = idx(2) + 1

            ! account for x periodicity
            call periodic_index_shift(is(:, n))

        end subroutine trilinear_interp

end module interpolation
