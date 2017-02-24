!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
program main
  use global_variables
  implicit none

  call initialize_mpi
  call set_model_parameters


  select case(calc_mode)
  case('MTEF')
    call allocation_of_global_arrays
    call multi_traject_Ehrenfest
  case('GQME_K')
    call allocation_of_global_arrays
    call GQME_kernel
    call GQME_dynamics
  case('GQME_T')
    call allocation_of_global_arrays
    call GQME_kernel_read
    call GQME_dynamics
  case('PBME')
    call PBME_allocation
    call PBME_dynamics
  case default
    call err_finalize('Invalid calc_mode')
  end select


  call MPI_finalize(ierr)

end program main
