!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GQME_dynamics
  use global_variables
  implicit none
  integer :: a1,a2
  integer :: it,it2
  complex(8),allocatable :: zrho_t(:,:),zrho_in(:,:),zrho_out(:,:)
  real(8) :: Ekin_s,norm_s
  if(myrank /= 0)return
  open(20,file="GQME_Et.out")

  call set_initial_conditions_elec
  allocate(zrho_DM(Lsite,Lsite,-1:Nt),zrho_t(Lsite,Lsite))
  allocate(zrho_in(Lsite,Lsite),zrho_out(Lsite,Lsite))

  do a1 = 1,Lsite
     do a2 = 1,Lsite
        zrho_DM(a1,a2,0) = zC(a1)*conjg(zC(a2))
     end do
  end do

  zrho_t=-zI*(matmul(Hmat_kin,zrho_DM(:,:,0))-matmul(zrho_DM(:,:,0),Hmat_kin))
  zrho_DM(:,:,-1) = zrho_DM(:,:,0) - dt*zrho_t

  zrho_t = matmul(Hmat_kin,zrho_DM(:,:,0))
  Ekin_s=0d0; norm_s=0d0
  do a1 = 1, Lsite
     Ekin_s = Ekin_s + zrho_t(a1,a1)
     norm_s = norm_s + zrho_DM(a1,a1,0)
  end do
  write(20,"(999e26.16e3)")0d0,norm_s,Ekin_s

  do it = 0,Nt-1

     zrho_DM(:,:,it+1) = zrho_DM(:,:,it-1)
     zrho_t=-zI*(matmul(Hmat_kin,zrho_DM(:,:,it))-matmul(zrho_DM(:,:,it),Hmat_kin))
     if(it /= 0)then
        call kernel_density_matrix(zK_full(:,:,:,:,0),zrho_DM(:,:,it),zrho_out,Lsite)
        zrho_t = zrho_t - 0.5d0*dt*zrho_out
        call kernel_density_matrix(zK_full(:,:,:,:,it),zrho_DM(:,:,0),zrho_out,Lsite)
        zrho_t = zrho_t - 0.5d0*dt*zrho_out
        do it2 = 1,it-1
           call kernel_density_matrix(zK_full(:,:,:,:,it2),zrho_DM(:,:,it-it2),zrho_out,Lsite)
           zrho_t = zrho_t - dt*zrho_out
        end do
     end if
     zrho_DM(:,:,it+1) = zrho_DM(:,:,it+1) + 2d0*dt*zrho_t

     zrho_t = matmul(Hmat_kin,zrho_DM(:,:,it+1))
     Ekin_s=0d0; norm_s=0d0
     do a1 = 1, Lsite
        Ekin_s = Ekin_s + zrho_t(a1,a1)
        norm_s = norm_s + zrho_DM(a1,a1,it+1)
     end do
     write(20,"(999e26.16e3)")dt*(it+1),norm_s,Ekin_s


  end do

  close(20)

end subroutine GQME_dynamics
!--------------------------------------------------------------
! calculate zK3 = zK1*zK2
subroutine kernel_density_matrix(zK1,zrho_in,zrho_out,Lsite)
  implicit none
  integer,intent(in) :: Lsite
  complex(8),intent(in) :: zK1(Lsite,Lsite,1,Lsite)
  complex(8),intent(in) :: zrho_in(Lsite,Lsite)
  complex(8),intent(out) :: zrho_out(Lsite,Lsite)
  complex(8) :: zs
  integer :: a1,a2,b1,b2,a1t,a2t,b1t,b2t

  do a1 = 1,Lsite
  do a2 = 1,Lsite

      zs = 0d0
      do b1 = 1,Lsite
        b1t=1
        a1t = mod((a1-b1)+ 2*Lsite,Lsite) + 1
        a2t = mod((a2-b1)+ 2*Lsite,Lsite) + 1
        do b2 = 1,Lsite
          b2t = mod((b2-b1)+ 2*Lsite,Lsite) + 1
          zs = zs + zK1(a1t,a2t,b1t,b2t)*zrho_in(b1,b2)
        end do
      end do

      zrho_out(a1,a2) = zs

  end do
  end do


end subroutine kernel_density_matrix
