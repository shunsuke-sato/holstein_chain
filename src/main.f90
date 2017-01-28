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

  call multi_traject_Ehrenfest


  call MPI_finalize(ierr)

end program main
