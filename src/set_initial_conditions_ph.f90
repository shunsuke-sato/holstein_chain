!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_initial_conditions_ph
  use global_variables
  implicit none

  call set_thermal_ph_dist 

  X_HO(:) = X_HO_ini(:)
  V_HO(:) = V_HO_ini(:)

end subroutine set_initial_conditions_ph
