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
  call refine_GQME_kernel_K1K3  ! Inversion symmetry, Unitarity, ... etc.
  call evaluate_GQME_kernel_full

  if(myrank == 0)then
    open(nfile_full_kernel,file=trim(file_full_kernel),form='unformatted')
    write(nfile_full_kernel)zK_full
    close(nfile_full_kernel)
  end if

end subroutine GQME_kernel
