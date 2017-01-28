!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine dt_evolve_elec
  use global_variables
  implicit none
  integer :: itraj,i
  complex(8) :: zCt_t(Lsite)
!LAPACK ==
  integer :: lwork,Nmat
  real(8) :: work(2*Lsite**2-2),diag(Lsite),off_diag(Lsite-1)
  real(8) :: eig(Lsite), eig_vec(Lsite,Lsite)
  integer :: info
!LAPACK ==

  
  diag=-gamma*sqrt(2d0*mass*omega0)*X_HO
  off_diag = - t0
  call dstev('V',Lsite,diag,off_diag,eig_vec,Lsite,work,info)

  zCt_t(:) = matmul(transpose(eig_vec(:,:)),zC(:))
  zC = 0d0
  do i = 1,Lsite
    zC(:)=zC(:) + exp(-zI*0.5d0*dt*diag(i))*zCt_t(i)*eig_vec(:,i)
  end do


  diag=-gamma*sqrt(2d0*mass*omega0)*X_HO_new
  off_diag = - t0
  call dstev('V',Lsite,diag,off_diag,eig_vec,Lsite,work,info)

  zCt_t(:) = matmul(transpose(eig_vec(:,:)),zC(:))
  zC = 0d0
  do i = 1,Lsite
    zC(:)=zC(:) + exp(-zI*0.5d0*dt*diag(i))*zCt_t(i)*eig_vec(:,i)
  end do


end subroutine dt_evolve_elec
