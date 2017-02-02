!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GQME_kernel
  use global_variables
  implicit none
  complex(8) :: zK_tmp(Lsite,Lsite),zK_tmp_l(Lsite,Lsite)
  real(8) :: phi
  integer :: itraj,it,jsite

  call allocate_GQME_kernel

! diagonal (1,1)
  do itraj = 1,Ntraj
    if(mod(itraj,Nprocs) /= myrank)cycle

    zC = 0d0; zC(1) = 1d0
    call set_initial_conditions_ph(itraj)
    call calc_force_HO
    X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

      call dt_evolve_elec
      X_HO_old = X_HO; X_HO = X_HO_new
      call calc_force_HO
! Kernel calc

    end do

  end do

! off-diagonal (1,j)
  do jsite = 2, Nsite

    do itraj = 1,Ntraj
      call random_number(phi); phi = 2d0 * pi
      if(mod(itraj,Nprocs) /= myrank)cycle

! positive summ
      zC = 0d0; zC(1) = sqrt(0.5d0); zC(jsite) = exp(zI*phi)*sqrt(0.5d0)
      call set_initial_conditions_ph(itraj)
      call calc_force_HO
      X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2
      
      do it = 0,Nt

         X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
         V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

         call dt_evolve_elec
         X_HO_old = X_HO; X_HO = X_HO_new
         call calc_force_HO
! Kernel calc
      end do

! negative summ
      zC = 0d0; zC(1) = sqrt(0.5d0); zC(jsite) = -exp(zI*phi)*sqrt(0.5d0)
      call set_initial_conditions_ph(itraj)
      call calc_force_HO
      X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2
      
      do it = 0,Nt

         X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
         V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

         call dt_evolve_elec
         X_HO_old = X_HO; X_HO = X_HO_new
         call calc_force_HO
! Kernel calc
      end do

    end do
  end do


end subroutine GQME_kernel
