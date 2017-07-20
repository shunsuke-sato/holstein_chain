!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine allocate_GQME_kernel
  use global_variables
  implicit none

  allocate(zK_full(Lsite,Lsite,1,Lsite,0:Nt),zK1(Lsite,Lsite,1,Lsite,0:Nt))
  allocate(zK2(Lsite,Lsite,1,Lsite,0:Nt),zK3(Lsite,Lsite,1,Lsite,0:Nt))


end subroutine Allocate_GQME_kernel
