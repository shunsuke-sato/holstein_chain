!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine multi_traject_Ehrenfest
  use global_variables
  implicit none
  integer :: itraj,it
  real(8) :: Ekin_t,Eph_t,Ecoup_t

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0

  do itraj = 1,Ntraj

    if(mod(itraj,Nprocs) /= myrank)cycle
    call set_initial_conditions_elec
    call set_initial_conditions_ph(itraj)
    call calc_energy(Ekin_t,Eph_t,Ecoup_t)
    Ekin_l(0)=Ekin_l(0)+Ekin_t
    Eph_l(0)=Eph_l(0)+Eph_t
    Ecoup_l(0)=Ecoup_l(0)+Ecoup_t

    do it = 0,Nt

      X_HO_new = X_HO_old + V_HO*2d0*dt + 0.5d0*F_HO/mass*(2d0*dt)**2
      call dt_evolve_elec
      call calc_force_HO
      V_HO_new = V_HO_old + (F_HO_new + 2d0*F_HO - F_HO_old)*dt

      X_HO = X_HO_new; V_HO = V_HO_new; F_HO = F_HO_new

      call calc_energy(Ekin_t,Eph_t,Ecoup_t)
      Ekin_l(it+1)=Ekin_l(it+1)+Ekin_t
      Eph_l(it+1)=Eph_l(it+1)+Eph_t
      Ecoup_l(it+1)=Ecoup_l(it+1)+Ecoup_t

    end do

  end do

  call MPI_ALLREDUCE(Ekin_l,Ekin,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Eph_l,Eph,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Ecoup_l,Ecoup,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  Ekin=Ekin/dble(Ntraj)
  Eph=Eph/dble(Ntraj)
  Ecoup=Ecoup/dble(Ntraj)

  if(myrank == 0)then
    open(21,file="energy_t.out")
    write(21,"(A)")"# tt, Etot, Ekin, Eph, Ecoup"
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,Ekin(it)+Eph(it)+Ecoup(it),Ekin(it),Eph(it),Ecoup(it)
    end do
    close(21)
  end if

end subroutine multi_traject_Ehrenfest
