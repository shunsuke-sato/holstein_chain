!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_dt_evolve
  use global_variables
  implicit none
  real(8) :: n_m(Lsite),n_n(Lsite),n(Lsite)
  real(8) :: Htot_t(Lsite,Lsite)
  integer :: i,j

  do i = 1,Lsite
    n_m(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    n_n(i) = 0.5d0*(x_n(i)**2 + p_n(i)**2)
    n(i) = 0.5d0*(n_m(i) + n_n(i)) - 1d0
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass


  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do

  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))
  p_n(:) = p_n(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_n(:))

  X_HO = X_HO + V_HO*dt

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do

  x_m(:) = x_m(:) + dt*matmul(Htot_t(:,:),p_m(:))
  x_n(:) = x_n(:) + dt*matmul(Htot_t(:,:),p_n(:))
  p_m(:) = p_m(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_m(:))
  p_n(:) = p_n(:) - 0.5d0*dt*matmul(Htot_t(:,:),x_n(:))

  do i = 1,Lsite
    n_m(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 )
    n_n(i) = 0.5d0*(x_n(i)**2 + p_n(i)**2 )
    n(i) = 0.5d0*(n_m(i) + n_n(i)) - 1d0
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

end subroutine FBTS_dt_evolve
