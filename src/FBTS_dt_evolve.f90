!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_dt_evolve_quantum
  use global_variables
  implicit none
  real(8) :: n_m(Lsite),n_n(Lsite),n(Lsite)
  real(8) :: Htot_t(Lsite,Lsite)
  complex(8) :: z_m(Lsite),z_n(Lsite)
  integer :: i,j
  complex(8) :: z_t(Lsite),zh_t(Lsite),zfact
  integer,parameter :: NTaylor=6
  integer :: iexp


  z_m = x_m + zI * p_m
  z_n = x_n + zI * p_n

  do i = 1,Lsite
    n_m(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    n_n(i) = 0.5d0*(x_n(i)**2 + p_n(i)**2)
!    n(i) = n_m(i) + n_n(i) ! with trace
    n(i) = n_m(i) + n_n(i) - 1.0d0 ! traceless
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2


  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do

  zfact = 1d0
  z_t = z_n
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_n = z_n + zfact*zh_t
    z_t = zh_t
  end do

  V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt + F_HO/mass*dt
  X_HO_old = X_HO; X_HO = X_HO_new

  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do

  zfact = 1d0
  z_t = z_n
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_n = z_n + zfact*zh_t
    z_t = zh_t
  end do


end subroutine FBTS_dt_evolve_quantum
!======================================================================================
subroutine FBTS_dt_evolve
  use global_variables
  implicit none
  real(8) :: n_m(Lsite),n_n(Lsite),n(Lsite)
  real(8) :: Htot_t(Lsite,Lsite)
  integer :: i,j

  do i = 1,Lsite
    n_m(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    n_n(i) = 0.5d0*(x_n(i)**2 + p_n(i)**2)
!    n(i) = n_m(i) + n_n(i) ! with trace
    n(i) = n_m(i) + n_n(i) - 1.0d0 ! traceless
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
!    n(i) = n_m(i) + n_n(i) ! with trace
    n(i) = n_m(i) + n_n(i) - 1.0d0 ! traceless
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

end subroutine FBTS_dt_evolve
!==========================================================================
subroutine FBTS_calc_ene(Etot0)
  use global_variables
  implicit none
  real(8) :: n_m(Lsite),n_n(Lsite),n(Lsite)
  real(8) :: Htot_t(Lsite,Lsite)
  integer :: i,j
  real(8) :: Etot0

  Etot0 = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)
  Etot0 = Etot0 - sum(  -gamma*sqrt(2d0*mass*omega0)*X_HO(:) ) ! traceless


  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do


  do i = 1,Lsite
    do j = 1,Lsite
      Etot0 = Etot0 + 0.5d0*Htot_t(i,j)*(x_m(i)*x_m(j) + p_m(i)*p_m(j) &
                                      +x_n(i)*x_n(j) + p_n(i)*p_n(j) )
    end do
  end do
 
end subroutine FBTS_calc_ene
