!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
program main
  use global_variables
  use CTEF_mod
  implicit none

  call initialize_mpi
  timer_start = MPI_Wtime()

  call initialize_random_number_generator
  call set_model_parameters


  select case(calc_mode)
  case('MTEF')
    call allocation_of_global_arrays
    call multi_traject_Ehrenfest
  case('MTEF_PAV')
    call allocation_of_global_arrays
    call multi_traject_Ehrenfest_phase_average
  case('GQME_K')
    call allocation_of_global_arrays
    call GQME_kernel
    call GQME_dynamics
  case('GQME_T')
    call allocation_of_global_arrays
    call GQME_kernel_read
    call GQME_dynamics
  case('PTEF')
    call PTEF_allocation
    call PTEF_dynamics
  case('CTEF')
    call CTEF
  case('CTEF_K')
    call CTEF_kernel
  case('PBME')
    call PBME_allocation
    call PBME_dynamics
  case('PBME_K')
    call PBME_allocation
    call PBME_kernel
    call GQME_dynamics
  case('FBTS')
    call FBTS_allocation
    call FBTS_dynamics
  case('FBTS_K')
    call FBTS_allocation
    call FBTS_kernel
    call GQME_dynamics
  case default
    call err_finalize('Invalid calc_mode')
  end select

  timer_end = MPI_Wtime()
  if(myrank == 0)write(*,"(A,2x,e16.6e3,2x,A)")"Total elapsed time =" &
    ,timer_end - timer_start,"sec."
  call MPI_finalize(ierr)

end program main
