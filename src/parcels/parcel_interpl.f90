! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : max_num_parcels, zero, one, two
    use parameters, only : nx, nz
    use options, only : parcel_info, interpl
    use parcel_container, only : parcel_container_type, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use ellipse
    use fields
    use taylorgreen, only : get_flow_velocity, &
                            get_flow_gradient

    implicit none

    private :: par2grid_elliptic,       &
               par2grid_non_elliptic,   &
               grid2par_elliptic,       &
               grid2par_non_elliptic,   &
               get_indices_and_weights


    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer :: is(ngp), js(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

    private :: is, js, weights

    contains

        subroutine par2grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)

            field = zero

            if (parcel_info%is_elliptic) then
                call par2grid_elliptic(parcels, attrib, field)
            else
                call par2grid_non_elliptic(parcels, attrib, field)
            endif

            ! apply free slip boundary condition
            field(0,  :, :) = two * field(0,  :, :)
            field(nz, :, :) = two * field(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            field(1,    :, :) = field(1,    :, :) + field(-1,   :, :)
            field(nz-1, :, :) = field(nz-1, :, :) + field(nz+1, :, :)

        end subroutine par2grid

        subroutine par2grid_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)
            integer                                    :: ncomp
            double precision                           :: points(2, 2)
            integer                                    :: n, p, c, l
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

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call get_indices_and_weights(points(p, :))

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do l = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            field(js(l), is(l), c) = field(js(l), is(l), c)         &
                                                   + 0.5d0 * weights(l) * attrib(n, c)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine par2grid_elliptic


        subroutine par2grid_non_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)
            integer                                    :: ncomp
            integer                                    :: n, c, l
            integer                                    :: the_shape(3)
            double precision                           :: pos(2)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                pos = parcels%position(n, :)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call get_indices_and_weights(pos)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        ! the weight is halved due to 2 points per ellipse
                        field(js(l), is(l), c) = field(js(l), is(l), c)   &
                                               + weights(l) * attrib(n, c)
                    enddo
                enddo
            enddo

        end subroutine par2grid_non_elliptic


        subroutine grid2par(position, volume, attrib, field, B, exact)
            double precision,           intent(in)  :: position(:, :)
            double precision,           intent(in)  :: volume(:, :)
            double precision,           intent(out) :: attrib(:, :)
            double precision,           intent(in)  :: field(-1:, 0:, :)
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
            double precision, intent(in)  :: field(-1:, 0:, :)
            logical, optional, intent(in) :: add
            character(*), optional, intent(in)      :: exact
            integer                       :: ncomp
            double precision              :: points(2, 2)
            integer                       :: n, p, c, l
            integer                       :: the_shape(3)


            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                   attrib(1:n_parcels, :) = zero
               endif
            else
               attrib(1:n_parcels, :) = zero
            endif

            ! put if statement here for computational efficiency
            if(present(exact)) then
               do n = 1, n_parcels

                  points = get_ellipse_points(position(n, :), volume(n, 1), B(n, :))

                  do p = 1, 2
                     call apply_periodic_bc(points(p, :))
                     if(exact=='velocity') then
                        attrib(n,:)=attrib(n,:) +0.5d0 * get_flow_velocity(points(p, :))
                     elseif(exact=='strain') then
                        attrib(n,:)=attrib(n,:) + 0.5d0 * get_flow_gradient(points(p, :))
                     else
                        print *, "Exact interpolation field passed not implemented"
                        stop
                     end if
                  end do
               end do
               return
            endif

            do n = 1, n_parcels

                points = get_ellipse_points(position(n, :), volume(n, 1), B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call get_indices_and_weights(points(p, :))

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do l = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            attrib(n, c) = attrib(n, c) &
                                         + 0.5d0 * weights(l) * field(js(l), is(l), c)
                        enddo
                    enddo
                enddo
            enddo

        end subroutine grid2par_elliptic


        subroutine grid2par_non_elliptic(position, attrib, field, add, exact)
            double precision, intent(in)  :: position(:, :)
            double precision, intent(out) :: attrib(:, :)
            double precision, intent(in)  :: field(-1:, 0:, :)
            logical, optional, intent(in) :: add
            character(*), optional, intent(in)      :: exact
            integer                       :: ncomp
            integer                       :: n, c, l
            integer                       :: the_shape(3)
            double precision              :: pos(2)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                   attrib(1:n_parcels, :) = zero
               endif
            else
               attrib(1:n_parcels, :) = zero
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

                pos = position(n, :)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call get_indices_and_weights(pos)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        attrib(n, c) = attrib(n, c) &
                                     + weights(l) * field(js(l), is(l), c)
                    enddo
                enddo
            enddo

        end subroutine grid2par_non_elliptic


        subroutine get_indices_and_weights(pos)
            double precision, intent(in)  :: pos(2)

            if (interpl == 'trilinear') then
                call trilinear(pos, is, js, weights)
            else if (interpl == 'exact') then ! only applies to par2grid
                call trilinear(pos, is, js, weights)
            else
                print *, "Unknown interpolation method '", interpl, "'."
                stop
            endif

        end subroutine get_indices_and_weights

        !
        ! tri-linear interpolation
        !
        subroutine trilinear(pos, ii, jj, ww)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ii(4), jj(4)
            double precision, intent(out) :: ww(4)
            double precision              :: xy(2)

            ! (i, j)
            call get_index(pos, ii(1), jj(1))
            call get_position(ii(1), jj(1), xy)
            ww(1) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j)
            ii(2) = ii(1) + 1
            jj(2) = jj(1)
            call get_position(ii(2), jj(2), xy)
            ww(2) = product(one - abs(pos - xy) * dxi)

            ! (i, j+1)
            ii(3) = ii(1)
            jj(3) = jj(1) + 1
            call get_position(ii(3), jj(3), xy)
            ww(3) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j+1)
            ii(4) = ii(2)
            jj(4) = jj(3)
            call get_position(ii(4), jj(4), xy)
            ww(4) = product(one - abs(pos - xy) * dxi)

            ! account for x periodicity
            call periodic_index_shift(ii)

        end subroutine trilinear

end module parcel_interpl
