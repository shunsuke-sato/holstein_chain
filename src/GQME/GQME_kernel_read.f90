!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GQME_kernel_read
  use global_variables
  implicit none

  call allocate_GQME_kernel

  if(myrank == 0)then
    open(nfile_full_kernel,file=trim(file_full_kernel),form='unformatted')
    read(nfile_full_kernel)zK_full,zK1,zK3
    close(nfile_full_kernel)
  end if

end subroutine GQME_kernel_read
