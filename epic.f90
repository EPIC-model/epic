program epic
    use model, only : pre_run, run, post_run
    implicit none

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

end program epic

