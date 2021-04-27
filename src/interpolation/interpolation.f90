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
    use interpl_methods
    use fields
    use taylorgreen, only : get_flow_velocity, &
                            get_flow_gradient

    implicit none

    private :: par2grid_elliptic,       &
               par2grid_non_elliptic,   &
               grid2par_elliptic,       &
               grid2par_non_elliptic,   &
               get_indices_and_weights


    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer ji(2, 4)

    ! interpolation weights
    double precision weight(4)

    private :: ji, weight

    contains

        subroutine par2grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)

            field = 0.0

            if (parcel_info%is_elliptic) then
                call par2grid_elliptic(parcels, attrib, field)
            else
                call par2grid_non_elliptic(parcels, attrib, field)
            endif

            ! apply free slip boundary condition
            field(0,  :, :) = 2.0 * field(0,  :, :)
            field(nz, :, :) = 2.0 * field(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            field(1,    :, :) = field(1,    :, :) + field(-1,   :, :)
            field(nz-1, :, :) = field(nz-1, :, :) + field(nz+1, :, :)

        end subroutine par2grid

        subroutine par2grid_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)
            integer                                    :: ncomp, ngp
            double precision                           :: points(2, 2)
            integer                                    :: n, p, c, i
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
                    call get_indices_and_weights(points(p, :), ngp)

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            field(ji(1, i), ji(2, i), c) = field(ji(1, i), ji(2, i), c)     &
                                                         + 0.5 * weight(i) * attrib(n, c)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine par2grid_elliptic


        subroutine par2grid_non_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(-1:, 0:, :)
            integer                                    :: ncomp, ngp
            integer                                    :: n, c, i
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
                call get_indices_and_weights(pos, ngp)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        ! the weight is halved due to 2 points per ellipse
                        field(ji(1, i), ji(2, i), c) = field(ji(1, i), ji(2, i), c) &
                                                     + weight(i) * attrib(n, c)
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
            integer                       :: ncomp, ngp
            double precision              :: points(2, 2)
            integer                       :: n, p, c, i
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

                points = get_ellipse_points(position(n, :), volume(n, 1), B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call get_indices_and_weights(points(p, :), ngp)

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            attrib(n, c) = attrib(n, c) &
                                         + 0.5 * weight(i) * field(ji(1, i), ji(2, i), c)
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

                pos = position(n, :)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call get_indices_and_weights(pos, ngp)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        attrib(n, c) = attrib(n, c) &
                                     + weight(i) * field(ji(1, i), ji(2, i), c)
                    enddo
                enddo
            enddo

        end subroutine grid2par_non_elliptic


        subroutine get_indices_and_weights(pos, ngp)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ngp


            if (interpl == 'trilinear') then
                call trilinear(pos, ji, weight, ngp)
            else if (interpl == 'exact') then ! only applies to par2grid
                call trilinear(pos, ji, weight, ngp)
            else
                print *, "Unknown interpolation method '", interpl, "'."
                stop
            endif

        end subroutine get_indices_and_weights

end module interpolation
