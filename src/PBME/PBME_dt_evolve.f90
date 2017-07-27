!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_dt_evolve_quantum
  use global_variables
  implicit none
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j
  complex(8) :: z_m(Lsite),z_t(Lsite),zh_t(Lsite),zfact
  integer,parameter :: NTaylor=6
  integer :: iexp

  z_m = x_m + zI * p_m

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do


  X_HO = X_HO + V_HO*dt

  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do

  x_m = real(z_m); p_m = aimag(z_m)

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do

  V_HO = V_HO +0.5d0*dt*F_HO/mass

end subroutine PBME_dt_evolve_quantum
!====================================================
subroutine PBME_C_dt_evolve_quantum
  use global_variables
  implicit none
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite),ntot,Eb1,Eb2
  integer :: i,j
  complex(8) :: z_m(Lsite),z_t(Lsite),zh_t(Lsite),zfact
  integer,parameter :: NTaylor=6
  integer :: iexp


  z_m = x_m + zI * p_m
  Eb1 = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i)
  end do
  ntot = sum(n)
  F_HO(:) = F_HO(:) - omega0**2*mass*X_HO(:) *ntot
  
  V_HO = V_HO +0.5d0*dt*F_HO/mass

  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do


  X_HO = X_HO + V_HO*dt*ntot

  zfact = 1d0
  z_t = z_m
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(z_t,zh_t)
    z_m = z_m + zfact*zh_t
    z_t = zh_t
  end do

  x_m = real(z_m); p_m = aimag(z_m)

  do i = 1,Lsite
    n(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2 - 1d0)
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i)
  end do
  ntot = sum(n)
  F_HO(:) = F_HO(:) - omega0**2*mass*X_HO(:) *ntot

  V_HO = V_HO +0.5d0*dt*F_HO/mass


  z_m = x_m + zI * p_m
  Eb2 = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)
  z_m = z_m*exp(-zI*dt*0.5d0*(Eb1+Eb2))
  x_m = real(z_m); p_m = aimag(z_m)

end subroutine PBME_C_dt_evolve_quantum
!====================================================
!====================================================
subroutine PBME_M_dt_evolve_quantum
  use global_variables
  implicit none
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite),ntot,Eb1,Eb2,Es
  real(8) :: V_HO_tmp(Lsite),Etot_t
  real(8) :: x2,weight,fact1,fact2,fact
  real(8) :: dt_t
  integer :: i,j,Nt_t,it_t
  complex(8) :: z_m(Lsite),z_t(Lsite),zh_t(Lsite),zfact
  integer,parameter :: NTaylor=6
  integer :: iexp

  x2 = sum(x_m**2 + p_m**2)
  weight = 2**(Lsite+1)*exp(-x2)
  fact1  = 2**(Lsite+1)*exp(-x2)
!  fact2  = 2**(Lsite+1)*exp(-x2)*abs(x2-dble(Lsite)*0.5d0)
!  fact = max(fact1,fact2)*10d0
  fact = fact1*20d0
  Nt_t = aint(fact)+1
  dt_t = dt/Nt_t

  do it_t = 1,Nt_t

    z_m = x_m + zI * p_m
    Eb1 = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)

    
    do i = 1,Lsite
      n(i) = weight*(x_m(i)**2 + p_m(i)**2 - 0.5d0)
      F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i)
    end do
    ntot = sum(n)

!    F_HO(:) = F_HO(:) - omega0**2*mass*X_HO(:) *ntot
    F_HO(:) = F_HO(:) - omega0**2*mass*X_HO(:)
    
    V_HO = V_HO +0.5d0*dt_t*F_HO/mass
    
    zfact = 1d0
    z_t = z_m
    do iexp = 1,NTaylor
      zfact = zfact*(-zI*dt_t*0.5d0)/iexp
      call hpsi(z_t,zh_t)
      zh_t = 2d0*weight*zh_t
      z_m = z_m + zfact*zh_t
      z_t = zh_t
    end do
    

    X_HO = X_HO + V_HO*dt_t
    
    zfact = 1d0
    z_t = z_m
    do iexp = 1,NTaylor
      zfact = zfact*(-zI*dt_t*0.5d0)/iexp
      call hpsi(z_t,zh_t)
      zh_t = 2d0*weight*zh_t
      z_m = z_m + zfact*zh_t
      z_t = zh_t
    end do
    
    x_m = real(z_m); p_m = aimag(z_m)
    
    do i = 1,Lsite
      n(i) = weight*(x_m(i)**2 + p_m(i)**2 - 0.5d0)
      F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i)
    end do
    F_HO(:) = F_HO(:) - omega0**2*mass*X_HO(:)! *ntot
    
    V_HO = V_HO +0.5d0*dt_t*F_HO/mass
    
!    z_m = x_m + zI * p_m
!    Eb2 = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)
!    z_m = z_m*exp(-zI*dt_t*0.5d0*(Eb1+Eb2))*exp(-zI*dt_t*(-2d0)*Etot_mod)

  end do
    
end subroutine PBME_M_dt_evolve_quantum
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
