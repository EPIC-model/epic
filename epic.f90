program epic
    use model, only : setup, run
    implicit none

    ! Create the model
    call setup

    ! Run the model
    call run

end program epic

