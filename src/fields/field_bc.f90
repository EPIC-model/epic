module field_bc
    use parameters, only : mesh
    implicit none

    private :: apply_free_slip,    &
               apply_periodic

    contains

        subroutine apply_field_bc(field)
            double precision, intent(inout) :: field(0:, 0:, :)
            integer                         :: i

            ! sum halo cell values to internal cells at the boundary
            field(1, :, :)            = field(1, :, :) + field(0, :, :)
            field(mesh%grid(1), :, :) = field(mesh%grid(1), :, :)      &
                                      + field(mesh%grid(1) + 1, :, :)

            field(:, 1, :)            = field(:, 1, :) + field(:, 0, :)
            field(:, mesh%grid(2), :) = field(:, mesh%grid(1), :)      &
                                      + field(:, mesh%grid(1) + 1, :)

            do i = 1, 2
                if (mesh%bc(i) == "free slip") then
                    call apply_free_slip(field, i)
                else if (mesh%bc(i) == "periodic") then
                    call apply_periodic(field, i)
                else
                    print *, "No boundary condition named '", mesh%bc(i), "'."
                    stop
                endif
            enddo
        end subroutine apply_field_bc


        subroutine apply_periodic(field, i)
            double precision, intent(inout) :: field(0:, 0:, :)
            integer,          intent(in)    :: i

            if (i == 1) then
                field(1, :, :)            = field(1, :, :)              &
                                          + field(mesh%grid(1), :, :)
                field(mesh%grid(1), :, :) = field(1, :, :)
            else if (i == 2) then
                field(:, 1, :)            = field(:, 1, :)              &
                                          + field(:, mesh%grid(2), :)
                field(:, mesh%grid(2), :) = field(:, 1, :)
            else
                print *, "Only up to 2 dimensions!"
            endif

        end subroutine apply_periodic


        subroutine apply_free_slip(field, i)
            double precision, intent(inout) :: field(0:, 0:, :)
            integer,          intent(in)    :: i

            if (i == 1) then
                field(1, :, :)            = 2.0 * field(1, :, :)
                field(mesh%grid(1), :, :) = 2.0 * field(mesh%grid(1), :, :)
            else if (i == 2) then
                field(:, 1, :)            = 2.0 * field(:, 1, :)
                field(:, mesh%grid(2), :) = 2.0 * field(:, mesh%grid(2), :)
            else
                print *, "Only up to 2 dimensions!"
            endif

        end subroutine apply_free_slip

end module field_bc
