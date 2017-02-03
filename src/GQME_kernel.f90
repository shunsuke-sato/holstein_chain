!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GQME_kernel
  use global_variables
  implicit none

  call allocate_GQME_kernel
  call evaluate_GQME_kernel_K1K3
  call evaluate_GQME_kernel_full

end subroutine GQME_kernel
