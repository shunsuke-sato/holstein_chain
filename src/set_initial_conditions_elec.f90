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
    zC(i) = exp(zI*2d0*pi*(i-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite))
  end do




end subroutine set_initial_conditions_elec
