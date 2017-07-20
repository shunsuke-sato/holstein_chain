!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine evaluate_GQME_kernel_full
  use global_variables
  implicit none
  integer,parameter :: Niter_scf = 20
  integer :: it,it2,iter_scf,i
  complex(8),allocatable :: zK2_tmp(:,:,:,:,:),zK_tmp(:,:,:,:),zK_sum(:,:,:,:)
  complex(8),allocatable :: zK2_l(:,:,:,:,:)
  integer :: mod_table(-Lsite:Lsite)

!  if(myrank /= 0)return

  allocate(zK2_tmp(Lsite,Lsite,1,Lsite,0:Nt),zK2_l(Lsite,Lsite,1,Lsite,0:Nt))
  allocate(zK_tmp(Lsite,Lsite,1,Lsite),zK_sum(Lsite,Lsite,1,Lsite))

  do i = -Lsite,Lsite
    mod_table(i) = mod(i+2*Lsite,Lsite)
  end do


! evaluate K2
  zK2 = zK3
  zK2_tmp = zK2
  if(myrank == 0)write(*,"(A,2x,2e16.6e3)")"K3/K1",sum(abs(zK3)**2)/sum(abs(zK1)**2)

  do iter_scf = 1,Niter_scf

    zK2_l = 0d0      
    do it = 0,Nt

      if(mod(it,Nprocs) /= myrank)cycle
      zK2_l(:,:,:,:,it) = zK3(:,:,:,:,it)
      if(it == 0)cycle
      call kernel_product(zK3(:,:,:,:,it),zK2_tmp(:,:,:,:,0),zK_tmp,Lsite,mod_table)
      zK_sum = 0.5d0*zK_tmp
      call kernel_product(zK3(:,:,:,:,0),zK2_tmp(:,:,:,:,it),zK_tmp,Lsite,mod_table)
      zK_sum = zK_sum + 0.5d0*zK_tmp
      do it2 = 1,it-1
        call kernel_product(zK3(:,:,:,:,it-it2),zK2_tmp(:,:,:,:,it2),zK_tmp,Lsite,mod_table)
        zK_sum = zK_sum + zK_tmp
      end do
      zK_sum = zK_sum*dt
      zK2_l(:,:,:,:,it) = zK2_l(:,:,:,:,it) + zI*zK_sum(:,:,:,:)
    end do

    call MPI_ALLREDUCE(zK2_l,zK2,(Nt+1)*Lsite**3, &
      MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(myrank == 0)write(*,"(A,2x,I5,2e16.6e3)")"K2 error",iter_scf &
      ,sum(abs(zK2-zK2_tmp)**2)/sum(abs(zK1)**2) &
         ,sum(abs(zK2)**2)/sum(abs(zK1)**2)
    zK2_tmp = zK2

  end do
  
! evaluate full-kernel
  zK2_l = 0d0
  do it = 0,Nt
      
    if(mod(it,Nprocs) /= myrank)cycle
    zK2_l(:,:,:,:,it) = zK1(:,:,:,:,it)
    if(it == 0)cycle
    call kernel_product(zK1(:,:,:,:,it),zK2(:,:,:,:,0),zK_tmp,Lsite,mod_table)
    zK_sum = 0.5d0*zK_tmp
    call kernel_product(zK1(:,:,:,:,0),zK2(:,:,:,:,it),zK_tmp,Lsite,mod_table)
    zK_sum = zK_sum + 0.5d0*zK_tmp
    do it2 = 1,it-1
      call kernel_product(zK1(:,:,:,:,it-it2),zK2(:,:,:,:,it2),zK_tmp,Lsite,mod_table)
      zK_sum = zK_sum + zK_tmp
    end do
    zK_sum = zK_sum*dt
    zK2_l(:,:,:,:,it) = zK2_l(:,:,:,:,it) + zI*zK_sum(:,:,:,:)
  end do

  call MPI_ALLREDUCE(zK2_l,zK_full,(Nt+1)*Lsite**3, &
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call refine_GQME_kernel_full

  if(myrank == 0)then
    open(20,file="k.out")
    do it = 0,Nt
      write(20,"(999e26.16e3)")dt*it,zK_full(1,2,1,2,it),conjg(zK_full(1,Lsite,1,Lsite,it))
    end do
    close(20)
    open(20,file="k1.out")
    do it = 0,Nt
      write(20,"(999e26.16e3)")dt*it,zK1(1,2,1,2,it),conjg(zK1(1,Lsite,1,Lsite,it))
    end do
    close(20)
    open(20,file="k3.out")
    do it = 0,Nt
      write(20,"(999e26.16e3)")dt*it,zK3(1,2,1,2,it),conjg(zK3(1,Lsite,1,Lsite,it))
    end do
    close(20)
  end if

end subroutine evaluate_GQME_kernel_full
!====================================================
! calculate zK3 = zK1*zK2
subroutine kernel_product(zK1,zK2,zK3,Lsite,mod_table)
  implicit none
  integer,intent(in) :: Lsite,mod_table(-Lsite:Lsite)
  complex(8),intent(in) :: zK1(Lsite,Lsite,1,Lsite),zK2(Lsite,Lsite,1,Lsite)
  complex(8),intent(out) :: zK3(Lsite,Lsite,1,Lsite)
  complex(8) :: zs
  integer :: a1,a2,b1,b2,c1,c2,a1t,a2t,c1t,c2t


!  zK3 = 0d0
!  do b1=1,1
    b1=1
    do b2=1,Lsite

      do a1 = 1,Lsite

        do a2 = 1,Lsite


          zs = 0d0
          do c2 = 1,Lsite
            do c1 = 1,Lsite
              c1t=1
              c2t = mod_table(c2-c1) + 1
              a1t = mod_table(a1-c1) + 1
              a2t = mod_table(a2-c1) + 1
              zs = zs + zK1(a1t,a2t,c1t,c2t)*zK2(c1,c2,b1,b2)
            end do
          end do

          zK3(a1,a2,b1,b2) =  zs
        end do
      end do
    end do


end subroutine kernel_product

!=============================================
subroutine kernel_product_org(zK1,zK2,zK3,Lsite,mod_table)
  implicit none
  integer,intent(in) :: Lsite,mod_table(-Lsite:Lsite)
  complex(8),intent(in) :: zK1(Lsite,Lsite,1,Lsite),zK2(Lsite,Lsite,1,Lsite)
  complex(8),intent(out) :: zK3(Lsite,Lsite,1,Lsite)
  complex(8) :: zs
  integer :: a1,a2,b1,b2,c1,c2,a1t,a2t,c1t,c2t

  do a1 = 1,Lsite
  do a2 = 1,Lsite
!  do b1=1,1
    b1=1
    do b2=1,Lsite

      zs = 0d0
      do c1 = 1,Lsite
        c1t=1
        a1t = mod_table(a1-c1) + 1
        a2t = mod_table(a2-c1) + 1
        do c2 = 1,Lsite
          c2t = mod_table(c2-c1) + 1
          zs = zs + zK1(a1t,a2t,c1t,c2t)*zK2(c1,c2,b1,b2)
        end do
      end do

      zK3(a1,a2,b1,b2) = zs

    end do
  end do
  end do


end subroutine kernel_product_org
