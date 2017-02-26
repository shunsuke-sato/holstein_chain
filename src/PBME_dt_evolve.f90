!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_dt_evolve
  use global_variables
  implicit none
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do

  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))

  X_HO = X_HO + V_HO*dt

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do

  x_m(:) = x_m(:) + dt*matmul(Htot_t(:,:),p_m(:))
  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

end subroutine PBME_dt_evolve
!====================================================
subroutine PBME_dt_evolve_traceless
  use global_variables
  implicit none
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  real(8) :: ss,xav
  integer :: i,j

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
    F_HO(i) = F_HO(i) - gamma*sqrt(2d0*mass*omega0)/Lsite
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

  Htot_t = Hmat_kin
  ss = 0d0
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
    ss = ss + Htot_t(i,i)
  end do
  ss = ss/Lsite
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - ss
  end do

  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))

  X_HO = X_HO + V_HO*dt


  Htot_t = Hmat_kin
  ss = 0d0
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
    ss = ss + Htot_t(i,i)
  end do
  ss = ss/Lsite
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - ss
  end do

  x_m(:) = x_m(:) + dt*matmul(Htot_t(:,:),p_m(:))
  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
    F_HO(i) = F_HO(i) - gamma*sqrt(2d0*mass*omega0)/Lsite
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

end subroutine PBME_dt_evolve_traceless
