!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine evaluate_GQME_kernel_full
  use global_variables
  implicit none
  integer,parameter :: Niter_scf = 25
  integer :: it,it2,iter_scf,i
  complex(8),allocatable :: zK2_tmp(:,:,:,:,:),zK_tmp(:,:,:,:),zK_sum(:,:,:,:)

  if(myrank /= 0)return

  allocate(zK2_tmp(Lsite,Lsite,1,Lsite,0:Nt))
  allocate(zK_tmp(Lsite,Lsite,1,Lsite),zK_sum(Lsite,Lsite,1,Lsite))



! evaluate K2
  zK2 = zK3
  zK2_tmp = zK2

  do iter_scf = 1,Niter_scf

    do it = 0,Nt
      
      zK2(:,:,:,:,it) = zK3(:,:,:,:,it)
      if(it == 0)cycle
      call kernel_product(zK3(:,:,:,:,it),zK2_tmp(:,:,:,:,0),zK_tmp,Lsite)
      zK_sum = 0.5d0*zK_tmp
      call kernel_product(zK3(:,:,:,:,0),zK2_tmp(:,:,:,:,it),zK_tmp,Lsite)
      zK_sum = zK_sum + 0.5d0*zK_tmp
      do it2 = 1,it-1
        call kernel_product(zK3(:,:,:,:,it-it2),zK2_tmp(:,:,:,:,it2),zK_tmp,Lsite)
        zK_sum = zK_sum + zK_tmp
      end do
      zK_sum = zK_sum*dt
      zK2(:,:,:,:,it) = zK2(:,:,:,:,it) + zI*zK_sum(:,:,:,:)
    end do

    write(*,"(A,2x,I5,e16.6e3)")"K2 error",iter_scf,sum(abs(zK2-zK2_tmp)**2)/sum(abs(zK2)**2)
    zK2_tmp = zK2

  end do
  
! evaluate full-kernel

    do it = 0,Nt
      
      zK_full(:,:,:,:,it) = zK1(:,:,:,:,it)
      if(it == 0)cycle
      call kernel_product(zK1(:,:,:,:,it),zK2(:,:,:,:,0),zK_tmp,Lsite)
      zK_sum = 0.5d0*zK_tmp
      call kernel_product(zK1(:,:,:,:,0),zK2(:,:,:,:,it),zK_tmp,Lsite)
      zK_sum = zK_sum + 0.5d0*zK_tmp
      do it2 = 1,it-1
        call kernel_product(zK1(:,:,:,:,it-it2),zK2(:,:,:,:,it2),zK_tmp,Lsite)
        zK_sum = zK_sum + zK_tmp
      end do
      zK_sum = zK_sum*dt
      zK_full(:,:,:,:,it) = zK_full(:,:,:,:,it) + zI*zK_sum(:,:,:,:)
    end do


end subroutine evaluate_GQME_kernel_full
!====================================================
! calculate zK3 = zK1*zK2
subroutine kernel_product(zK1,zK2,zK3,Lsite)
  implicit none
  integer,intent(in) :: Lsite
  complex(8),intent(in) :: zK1(Lsite,Lsite,1,Lsite),zK2(Lsite,Lsite,1,Lsite)
  complex(8),intent(out) :: zK3(Lsite,Lsite,1,Lsite)
  complex(8) :: zs
  integer :: a1,a2,b1,b2,c1,c2,c1t,c2t

  do a1 = 1,Lsite
  do a2 = 1,Lsite
!  do b1=1,1
    b1=1
    do b2=1,Lsite

      zs = 0d0
      do c1 = 1,Lsite
        c1t=1
        do c2 = 1,Lsite
          c2t = mod((c2-c1)+ 2*Lsite,Lsite) + 1
          zs = zs + zK1(a1,a2,c1t,c2t)*zK2(c1,c2,b1,b2)
        end do
      end do

      zK3(a1,a2,b1,b2) = zs

    end do
  end do
  end do


end subroutine kernel_product
