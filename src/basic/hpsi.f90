!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine hpsi(zCt_t,zhCt_t)
  use global_variables
  implicit none
  integer :: i
  complex(8),intent(in) :: zCt_t(Lsite)
  complex(8),intent(out) :: zhCt_t(Lsite)


  
  zhCt_t=-gamma*sqrt(2d0*mass*omega0)*X_HO*zCt_t

  zhCt_t(1) = zhCt_t(1) -t0*(zCt_t(2) + zCt_t(Lsite))
  do i = 2,Lsite-1
    zhCt_t(i) = zhCt_t(i) -t0*(zCt_t(i+1) + zCt_t(i-1))
  end do
  zhCt_t(Lsite) = zhCt_t(Lsite) -t0*(zCt_t(1) + zCt_t(Lsite-1))


end subroutine hpsi
!=======================================================================
subroutine hpsi_new(zCt_t,zhCt_t)
  use global_variables
  implicit none
  integer :: i
  complex(8),intent(in) :: zCt_t(Lsite)
  complex(8),intent(out) :: zhCt_t(Lsite)


  
  zhCt_t=-gamma*sqrt(2d0*mass*omega0)*X_HO_new*zCt_t

  zhCt_t(1) = zhCt_t(1) -t0*(zCt_t(2) + zCt_t(Lsite))
  do i = 2,Lsite-1
    zhCt_t(i) = zhCt_t(i) -t0*(zCt_t(i+1) + zCt_t(i-1))
  end do
  zhCt_t(Lsite) = zhCt_t(Lsite) -t0*(zCt_t(1) + zCt_t(Lsite-1))


end subroutine hpsi_new
