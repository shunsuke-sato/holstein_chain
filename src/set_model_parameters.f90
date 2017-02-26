!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_model_parameters
  use global_variables
  implicit none

  Lsite = 12
  t0 = 1d0
  omega0 = 1d0
  gamma = sqrt(0.4d0)
  mass = 1d0

  Tph = -1d0

  Ntraj = 1000

  dt = 0.01d0
!  Nt = aint(30d0/dt)+1
  Nt = aint(25d0/dt)+1
!  Nt = aint(10d0/dt)+1

!'MTEF', 'GQME_K'

!  calc_mode = 'MTEF'
!  calc_mode = 'GQME_K'
!  calc_mode = 'GQME_T'
  calc_mode = 'PBME'



end subroutine set_model_parameters
