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


    call set_initial_conditions
    if(mod(itraj,Nprocs) /= myrank)cycle

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      call calc_energy(Ekin_t,Eph_t,Ecoup_t)
      Ekin_l(it)=Ekin_l(it)+Ekin_t
      Eph_l(it)=Eph_l(it)+Eph_t
      Ecoup_l(it)=Ecoup_l(it)+Ecoup_t

      call dt_evolve_elec

      X_HO_old = X_HO; X_HO = X_HO_new
      call calc_force_HO

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
    do it = 0,Nt
      write(21,"(999e26.16e3)")dt*it,Ekin(it)+Eph(it)+Ecoup(it),Ekin(it),Eph(it),Ecoup(it)
    end do
    close(21)
  end if

end subroutine multi_traject_Ehrenfest
