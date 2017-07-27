!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine refine_GQME_kernel_K1K3
  use global_variables
  implicit none
  complex(8),allocatable :: zK1_tmp(:,:,:,:,:),zK3_tmp(:,:,:,:,:)
  integer :: a1,a2,b1,b2,a1i,a2i,b1i,b2i,iter

!  return
  allocate(zK1_tmp(Lsite,Lsite,1,Lsite,0:Nt),zK3_tmp(Lsite,Lsite,1,Lsite,0:Nt))

!! Inversion symmetry
!  do a1 = 1,Lsite
!    a1i = mod(-(a1-1)+2*Lsite,Lsite)+1
!    do a2 = 1,Lsite
!      a2i = mod(-(a2-1)+2*Lsite,Lsite)+1
!      do b2 = 1,Lsite
!        b2i = mod(-(b2-1)+2*Lsite,Lsite)+1
!    
!        zK1_tmp(a1,a2,1,b2,:) = 0.5d0*(zK1(a1,a2,1,b2,:) + zK1(a1i,a2i,1,b2i,:) )
!        zK3_tmp(a1,a2,1,b2,:) = 0.5d0*(zK3(a1,a2,1,b2,:) + zK3(a1i,a2i,1,b2i,:) )
!        
!      end do
!    end do
!  end do
!
!  zK1 = zK1_tmp
!  zK3 = zK3_tmp


! self-joint symmetry
  do a1 = 1,Lsite
    do a2 = 1,Lsite
      do b2 = 1,Lsite
        a1i = mod(a1-(b2-1) -1 + 2*Lsite,Lsite) +1
        a2i = mod(a2-(b2-1) -1 + 2*Lsite,Lsite) +1
        b1i = mod(1-(b2-1) -1 + 2*Lsite,Lsite) +1
    
        zK1_tmp(a1,a2,1,b2,:) = 0.5d0*(zK1(a1,a2,1,b2,:) + conjg(zK1(a2i,a1i,1,b1i,:)) )
        zK3_tmp(a1,a2,1,b2,:) = 0.5d0*(zK3(a1,a2,1,b2,:) - conjg(zK3(a2i,a1i,1,b1i,:)) )
        
      end do
    end do
  end do

  zK1 = zK1_tmp
  zK3 = zK3_tmp


end subroutine refine_GQME_kernel_K1K3
