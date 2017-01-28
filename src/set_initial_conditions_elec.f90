!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_initial_conditions_elec
  use global_variables
  implicit none
  integer :: i

  do i = 1,Lsite
    zC(i) = exp(zI*pi*(i-1))/sqrt(dble(Lsite))
  end do




end subroutine set_initial_conditions_elec
