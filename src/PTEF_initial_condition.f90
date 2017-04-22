!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PTEF_initial_condition
  use global_variables
  implicit none
  integer :: isite,jsite
  real(8) :: xx,pp,ss

  call set_initial_conditions_elec
  zCp = zC
  call set_initial_conditions_elec

  do isite = 1,Lsite
    call correlated_gaussian_random_number(xx,pp)
    X_HO(isite) = sqrt(2d0/(mass*omega0))*xx
    Xp_HO(isite) = sqrt(2d0/(mass*omega0))*pp
    call correlated_gaussian_random_number(xx,pp)
    V_HO(isite) = sqrt(2d0*mass*omega0)*xx
    Vp_HO(isite) = sqrt(2d0*mass*omega0)*pp
  end do

  call calc_force_HO_pair
  X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2
  Xp_HO_old = Xp_HO - Vp_HO*dt +0.5d0*Fp_HO/mass*dt**2

end subroutine PTEF_initial_condition
