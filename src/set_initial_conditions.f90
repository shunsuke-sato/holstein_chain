!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_initial_conditions(itraj)
  use global_variables
  implicit none
  integer, intent(in) :: itraj

  call set_initial_conditions_elec
  call set_initial_conditions_ph(itraj)

  call calc_force_HO
  X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

end subroutine set_initial_conditions
