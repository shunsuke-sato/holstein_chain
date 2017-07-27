!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine dt_evolve_elec
  use global_variables
  implicit none
  integer,parameter :: NTaylor=6
  integer :: iexp
  complex(8) :: zCt_t(Lsite),zhCt_t(Lsite),zfact


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



end subroutine dt_evolve_elec
