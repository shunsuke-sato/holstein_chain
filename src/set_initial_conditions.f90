!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_initial_conditions
  use global_variables
  implicit none

  call set_initial_conditions_elec
  call set_initial_conditions_ph

  call calc_force_HO
  F_HO = F_HO_new;F_HO_old = F_HO_new
  V_HO_old = V_HO -F_HO/mass*dt

end subroutine set_initial_conditions
