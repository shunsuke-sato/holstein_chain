!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine calc_force_HO
  use global_variables
  implicit none
  integer :: i

  do i = 1,Lsite
    F_HO_new(i) = gamma*abs(zC(i))**2*sqrt(2d0*mass*omega0)
  end do
  

end subroutine calc_force_HO
