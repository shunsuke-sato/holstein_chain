!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine initialize_mpi
  use global_variables
  use parallel
  implicit none

  call init_parallel

  Nprocs = comm_nproc_global
  Myrank = comm_id_global

end subroutine initialize_mpi
