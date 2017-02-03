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
  call allocation_of_global_arrays
  call set_thermal_ph_dist

  select case(calc_mode)
  case('MTEF')
    call multi_traject_Ehrenfest
  case('GQME_K')
    call GQME_kernel
    call GQME_dynamics
  case default
    call err_finalize('Invalid calc_mode')
  end select


  call MPI_finalize(ierr)

end program main
