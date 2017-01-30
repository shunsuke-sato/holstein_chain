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

  Ntraj = 2

  Nt = 8000
  dt = 0.01d0

end subroutine set_model_parameters
