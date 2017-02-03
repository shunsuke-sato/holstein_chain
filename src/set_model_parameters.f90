!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_model_parameters
  use global_variables
  implicit none

  Lsite = 8
  t0 = 1d0
  omega0 = 1d0
  gamma = 1d0
  mass = 1d0

  Tph = -1d0

  Ntraj = 400

  dt = 0.08d0
!  Nt = aint(30d0/dt)+1
  Nt = aint(5d0/dt)+1

!'MTEF', 'GQME_K'
  calc_mode = 'GQME_K'

end subroutine set_model_parameters
