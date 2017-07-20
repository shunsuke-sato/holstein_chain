!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine dt_evolve_elec_pair
  use global_variables
  implicit none
  integer,parameter :: NTaylor=6
  integer :: iexp
  complex(8) :: zCt_t(Lsite),zhCt_t(Lsite),zfact
  integer,parameter :: iswitch = 1


  zfact = 1d0
  zCt_t = zC
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi(zCt_t,zhCt_t)
    zC = zC + zfact*zhCt_t
    zCt_t = zhCt_t
  end do


  zfact = 1d0
  zCt_t = zC
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi_new(zCt_t,zhCt_t)
    zC = zC + zfact*zhCt_t
    zCt_t = zhCt_t
  end do



  zfact = 1d0
  zCt_t = zCp
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi_pair(zCt_t,zhCt_t)
    zCp = zCp + zfact*zhCt_t
    zCt_t = zhCt_t
  end do


  zfact = 1d0
  zCt_t = zCp
  do iexp = 1,NTaylor
    zfact = zfact*(-zI*dt*0.5d0)/iexp
    call hpsi_new_pair(zCt_t,zhCt_t)
    zCp = zCp + zfact*zhCt_t
    zCt_t = zhCt_t
  end do



end subroutine dt_evolve_elec_pair

!=======================================================================
subroutine hpsi_pair(zCt_t,zhCt_t)
  use global_variables
  implicit none
  integer :: i
  complex(8),intent(in) :: zCt_t(Lsite)
  complex(8),intent(out) :: zhCt_t(Lsite)


  
  zhCt_t=-gamma*sqrt(2d0*mass*omega0)*Xp_HO*zCt_t
!  zhCt_t=zhCt_t + 1d0*zCt_t ! for test

  zhCt_t(1) = zhCt_t(1) -t0*(zCt_t(2) + zCt_t(Lsite))
  do i = 2,Lsite-1
    zhCt_t(i) = zhCt_t(i) -t0*(zCt_t(i+1) + zCt_t(i-1))
  end do
  zhCt_t(Lsite) = zhCt_t(Lsite) -t0*(zCt_t(1) + zCt_t(Lsite-1))


end subroutine hpsi_pair
!=======================================================================
subroutine hpsi_new_pair(zCt_t,zhCt_t)
  use global_variables
  implicit none
  integer :: i
  complex(8),intent(in) :: zCt_t(Lsite)
  complex(8),intent(out) :: zhCt_t(Lsite)


  
  zhCt_t=-gamma*sqrt(2d0*mass*omega0)*Xp_HO_new*zCt_t
!  zhCt_t=zhCt_t + 1d0*zCt_t ! for test

  zhCt_t(1) = zhCt_t(1) -t0*(zCt_t(2) + zCt_t(Lsite))
  do i = 2,Lsite-1
    zhCt_t(i) = zhCt_t(i) -t0*(zCt_t(i+1) + zCt_t(i-1))
  end do
  zhCt_t(Lsite) = zhCt_t(Lsite) -t0*(zCt_t(1) + zCt_t(Lsite-1))


end subroutine hpsi_new_pair
