!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_initial_conditions_ph(itraj)
  use global_variables
  implicit none
  integer, intent(in) :: itraj

  X_HO(:) = X_HO_ini(:,itraj)
  V_HO(:) = V_HO_ini(:,itraj)

end subroutine set_initial_conditions_ph
