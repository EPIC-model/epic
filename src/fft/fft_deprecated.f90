!     use datatypes, only : DEFAULT_PRECISION
    !use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, global_grid_type
    !use state_mod, only : model_state_type
    !use mpi, only : MPI_DOUBLE_COMPLEX, MPI_INT, MPI_COMM_SELF
    !use omp_lib
    !use ffte_mod, only: ffte_r2c, ffte_c2r, ffte_init, ffte_finalise

!> Performs a forward 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
!! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
!! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and
!! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
!! @param current_state The current model state
!! @param source_data The source real data to in the time domain
!! @param target_data Frequency domain real representation of the time domain source which is allocated here
subroutine perform_forward_3dfft(current_state, source_data, target_data)
type(model_state_type), target, intent(inout) :: current_state
double precision, dimension(:,:,:), intent(inout) :: source_data
double precision, dimension(:,:,:), intent(out) :: target_data

call transpose_and_forward_fft_in_y(current_state, source_data, buffer1, real_buffer1)

call transpose_and_forward_fft_in_x(current_state, real_buffer1, buffer2, real_buffer2)


call transpose_to_pencil(y_from_x_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
        real_buffer2, real_buffer3)
call transpose_to_pencil(z_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
    real_buffer3, target_data)


end subroutine perform_forward_3dfft

!> Performs a backwards 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
!! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
!! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and
!! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
!! @param current_state The current model state
!! @param source_data The source real data to in the frequency domain
!! @param target_data Time domain complex representation of the frequency domain source
subroutine perform_backwards_3dfft(current_state, source_data, target_data)
type(model_state_type), target, intent(inout) :: current_state
double precision, dimension(:,:,:), intent(in) :: source_data
double precision, dimension(:,:,:), intent(out) :: target_data

call transpose_to_pencil(y_from_z_2_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
    source_data, real_buffer3)
call transpose_to_pencil(x_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
    real_buffer3, real_buffer2)

call transpose_and_backward_fft_in_x(current_state, real_buffer2, buffer2, real_buffer1)
call transpose_and_backward_fft_in_y(current_state, real_buffer1, buffer1, target_data)

end subroutine perform_backwards_3dfft


!> Performs the transposition and forward FFT in the y dimension then converts back to real numbers. The Y size is
!! (n/2+1)*2 due to the complex to real transformation after the FFT.
!! @param current_state The current model state
!! @param source_data Input buffer, Z pencil oriented z,y,x
!! @param buffer Complex buffer which the FFT writes into
!! @param real_buffer Output buffer, Y pencil, oriented y,x,z
subroutine transpose_and_forward_fft_in_y(current_state, source_data, buffer, real_buffer)
type(model_state_type), target, intent(inout) :: current_state
double precision, dimension(:,:,:), intent(inout) :: source_data
double precision, dimension(:,:,:),  intent(out) :: real_buffer
complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer

! Transpose globally from Z pencil to Y pencil
call transpose_to_pencil(y_from_z_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
    source_data, fft_in_y_buffer)

call perform_r2c_fft(fft_in_y_buffer, buffer, y_from_z_transposition%my_pencil_size(Y_INDEX), &
        y_from_z_transposition%my_pencil_size(X_INDEX) * y_from_z_transposition%my_pencil_size(Z_INDEX), 1)
call convert_complex_to_real(buffer, real_buffer)
end subroutine transpose_and_forward_fft_in_y

!> Performs the backwards FFT in X and then transposes to Y pencil. The FFT requires complex numbers which are converted to real,
!! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
!! size is n+2, complex number size is n/2+1 and we get n reals out.
!! @param current_state The current model state
!! @param source_data Input buffer, X pencil oriented x,z,y
!! @param buffer Complex buffer which is fed into the FFT
!! @param real_buffer Output buffer, Y pencil, oriented y,x,z
subroutine transpose_and_backward_fft_in_x(current_state, source_data, buffer, real_buffer)
type(model_state_type), target, intent(inout) :: current_state
double precision, dimension(:,:,:), intent(inout) :: source_data
double precision, dimension(:,:,:),  intent(out) :: real_buffer
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

call convert_real_to_complex(source_data, buffer)
call perform_c2r_fft(buffer, fft_in_x_buffer, x_from_y_2_transposition%my_pencil_size(X_INDEX)-2, &
        x_from_y_2_transposition%my_pencil_size(Y_INDEX) * x_from_y_2_transposition%my_pencil_size(Z_INDEX), 2)

! Transpose globally from X pencil to Y pencil
call transpose_to_pencil(y_from_x_2_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
    fft_in_x_buffer, real_buffer)
end subroutine transpose_and_backward_fft_in_x

!> Performs the transposition and forward FFT in the x dimension. After the FFT the complex space is converted back into
!! real numbers. The X size is (n/2+1)*2 due to this transformation.
!! @param current_state The current model state
!! @param buffer1 Input buffer, Y pencil after the Y dimension FFT oriented y,x,z
!! @param buffer Complex buffer which results from the FFT
!! @param buffer2 Output buffer, X pencil after this X FFT, oriented x,z,y
subroutine transpose_and_forward_fft_in_x(current_state, source_data, buffer, real_buffer)
type(model_state_type), target, intent(inout) :: current_state
complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer
double precision, dimension(:,:,:), intent(inout) :: source_data, real_buffer

! Go from global Y pencil to global X pencil
call transpose_to_pencil(x_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
    source_data, fft_in_x_buffer)

call perform_r2c_fft(fft_in_x_buffer, buffer, x_from_y_transposition%my_pencil_size(X_INDEX), &
        x_from_y_transposition%my_pencil_size(Y_INDEX) * x_from_y_transposition%my_pencil_size(Z_INDEX), 3)

call convert_complex_to_real(buffer, real_buffer)
end subroutine transpose_and_forward_fft_in_x


!> Performs the backwards FFT in Y and then transposes to Z pencil. The FFT requires complex numbers which are converted to real,
!! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
!! size is n+2, complex number size is n/2+1 and we get n reals out.
!! @param current_state The current model state
!! @param source_data Input buffer, Y pencil oriented y,x,z
!! @param buffer Complex buffer which is fed into the FFT
!! @param real_buffer Output buffer, Z pencil, oriented z,y,x
subroutine transpose_and_backward_fft_in_y(current_state, source_data, buffer, real_buffer)
type(model_state_type), target, intent(inout) :: current_state
double precision, dimension(:,:,:), intent(inout) :: source_data
double precision, dimension(:,:,:),  intent(out) :: real_buffer
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

call convert_real_to_complex(source_data, buffer)

call perform_c2r_fft(buffer, fft_in_y_buffer,  y_from_x_2_transposition%my_pencil_size(Y_INDEX)-2, &
        y_from_x_2_transposition%my_pencil_size(X_INDEX) * y_from_x_2_transposition%my_pencil_size(Z_INDEX), 4)

! Go from global Y pencil to global Z pencil
call transpose_to_pencil(z_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
    fft_in_y_buffer, real_buffer)
end subroutine transpose_and_backward_fft_in_y



!> Actually performs a forward real to complex FFT
!! @param source_data Source (real) data in the time domain
!! @param transformed_data Resulting complex data in the frequency domain
!! @param row_size Number of elements for each FFT
!! @param num_rows The number of FFTs to perform on the next data elements in the source_data
!! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
subroutine perform_r2c_fft(source_data, transformed_data, row_size, num_rows, plan_id)
double precision, dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
integer, intent(in) :: row_size, num_rows, plan_id
integer :: i, j

! if (ffte) then !use FFTE for the FFTs

    !$OMP SINGLE
    call ffte_init(row_size)
    !$OMP END SINGLE

    !$OMP DO private(j)
    do i=1,size(source_data,3)
    do j=1,size(source_data,2)
        call ffte_r2c(source_data(:,j,i),transformed_data(:,j,i),row_size)
    enddo
    enddo
    !$OMP END DO

    !make sure all the threads have completed the above do loops before finalising
    !$OMP BARRIER

    !$OMP SINGLE
    call ffte_finalise()

    !$OMP END SINGLE


end subroutine perform_r2c_fft

!> Performs the complex to real (backwards) FFT
!! @param source_data Source (complex) data in the frequency domain
!! @param transformed_data Resulting real data in the time domain
!! @param row_size Number of elements for each FFT
!! @param num_rows The number of FFTs to perform on the next data elements in the source_data
!! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
subroutine perform_c2r_fft(source_data, transformed_data, row_size, num_rows, plan_id)
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
double precision, dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
integer, intent(in) :: row_size, num_rows, plan_id
integer :: i,j

!$OMP SINGLE
call ffte_init(row_size)
!$OMP END SINGLE

!$OMP DO private(j)
do i=1,size(source_data,3)
    do j=1,size(source_data,2)
    call ffte_c2r(source_data(:,j,i),transformed_data(:,j,i),row_size)
    enddo
enddo
!$OMP END DO

!make sure all the threads have completed the above do loops before finalising
!$OMP BARRIER

!$OMP SINGLE
call ffte_finalise()

!$OMP END SINGLE

end subroutine perform_c2r_fft


!> Converts complex representation to its real data counterpart and is called after each forward FFT.
!! After a r2c FFT, there are n/2+1 complex numbers - which means that there will be more real numbers in Fourier space
!! than are provided into the forward FFT call (due to the extra +1). Note that the real size n will always be complex size * 2
!! This always unpacks the complex dimension in the first dimension
!! @param complex_data Complex data in Z,Y,X orientation to be unpacked into its real representation
!! @param real_data The real representation is written into here
subroutine convert_complex_to_real(complex_data, real_data)
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), intent(in) :: complex_data
double precision, dimension(:,:,:), intent(out) :: real_data

integer :: i, j, k

!$OMP DO
do i=1,size(real_data,3)
    do j=1,size(real_data,2)
    do k=1,size(real_data,1),2
        real_data(k,j,i)=real(real(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
        real_data(k+1,j,i)=real(aimag(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
    end do
    end do
end do
!$OMP END DO


end subroutine convert_complex_to_real

!> Converts reals into their complex representation, this is called for backwards FFTs as we need to feed in complex numbers
!! to force FFTE to do a backwards. It is a relatively simple transformation, as n goes into n/2 complex numbers and as this
!! is the result of the `convert_complex_to_real` procedure, n always divides evenly.
!! This is always applied to the first dimension of the real data
!! @param real_data The source real data to pack into the complex data, it is oriented Z,Y,X
!! @param complex_data Target complex data which the real data is packaged into
subroutine convert_real_to_complex(real_data, complex_data)
double precision, dimension(:,:,:), intent(in) :: real_data
complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: complex_data

integer :: i, j, k

!$OMP WORKSHARE
complex_data(:,:,:)=cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
!$OMP END WORKSHARE

!$OMP DO
do i=1,size(real_data,3)
    do j=1,size(real_data,2)
    do k=1,size(real_data,1),2
        complex_data((k+1)/2,j,i)=cmplx(real_data(k,j,i), real_data(k+1,j,i), kind=C_DOUBLE_COMPLEX)
    end do
    end do
end do
!$OMP END DO

end subroutine convert_real_to_complex

!> Determines my global start coordinate in Fourier space.
!! This is required for cos y and cos x calculation which is fed into the tridiagonal solver. After the forward FFTs,
!! each process has ((n/2+1)/p+r) * 2 elements, where p is the number of processes and r
!! is the uneven process remainder (1 or 0 depending on p). Therefore some processes will have t elements, and some t-2 elements
!! to feed into the solver
!! @param current_state The current model state
!! @param dimension The dimension that we are calculating this for (Y or X)
!! @returns My global start in Fourier space
integer function deduce_my_global_start(current_state, dimension)
type(model_state_type), intent(inout) :: current_state
integer, intent(in) :: dimension

integer complex_size, distributed_size, remainder, larger_nums, smaller_nums

complex_size=(current_state%global_grid%size(dimension)/2+1)*2
distributed_size=complex_size / current_state%parallel%dim_sizes(dimension)
remainder=complex_size - distributed_size * current_state%parallel%dim_sizes(dimension)
larger_nums=min(remainder, current_state%parallel%my_coords(dimension))
smaller_nums=current_state%parallel%my_coords(dimension)-remainder
deduce_my_global_start=((distributed_size+1)*larger_nums + merge(distributed_size*smaller_nums, 0, smaller_nums .gt. 0)) + 1
end function deduce_my_global_start
