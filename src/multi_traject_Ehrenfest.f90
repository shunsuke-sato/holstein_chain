!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine multi_traject_Ehrenfest
  use global_variables
  implicit none
  integer :: itraj,it

  do itraj = 1,Ntraj

    if(mod(itraj,Nprocs) /= myrank)cycle
    call set_initial_conditions_elec
    call set_initial_conditions_ph(itraj)

    do it = 0,Nt

      X_HO_new = X_HO_old + V_HO*2d0*dt + 0.5d0*F_HO/mass*(2d0*dt)**2
      call dt_evolve_elec
      call calc_force_HO
      V_HO_new = V_HO_old + (F_HO_new + 2d0*F_HO - F_HO_old)*dt

    end do

  end do


end subroutine multi_traject_Ehrenfest
