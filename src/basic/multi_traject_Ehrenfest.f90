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
!-----------------------------------------------------------------------------------------
subroutine multi_traject_Ehrenfest_phase_average
  use global_variables
  implicit none
  integer :: itraj,it, i,j
  real(8) :: Ekin_t,Eph_t,Ecoup_t,phi
  complex(8) :: zrho_dm_fact

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0

  do itraj = 1,Ntraj


    zC = 0d0
    i = 1; j = mod(itraj,Lsite)+1
    if(i == j)then
      zC(i) = 1d0
    else
      call random_number(phi); phi = 2d0*phi*pi
      zC(i) = sqrt(0.5d0)
      zC(j) = exp(zI*phi)*sqrt(0.5d0)
    end if
    zrho_dm_fact = exp(zI*2d0*pi*(j-i)*dble(Lsite/2)/dble(Lsite))/dble(Lsite)*dble(Lsite**2)
    if(i /= j)zrho_dm_fact = zrho_dm_fact*exp(-zI*phi)*2d0
    call set_initial_conditions_ph
    call calc_force_HO
    X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2
    if(mod(itraj,Nprocs) /= myrank)cycle

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      call calc_energy(Ekin_t,Eph_t,Ecoup_t)
      Ekin_l(it)=Ekin_l(it)+Ekin_t*zrho_dm_fact
      Eph_l(it)=Eph_l(it)+Eph_t*zrho_dm_fact
      Ecoup_l(it)=Ecoup_l(it)+Ecoup_t*zrho_dm_fact

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
    open(21,file="MTEF_phase_ave.out")
    write(21,"(A)")"# tt, Ekin, Eph, Ecoup, Etot"
    do it = 0,Nt
      write(21,"(999e26.16e3)")dt*it,Ekin(it),Eph(it),Ecoup(it),Ekin(it)+Eph(it)+Ecoup(it)
    end do
    close(21)
  end if

end subroutine multi_traject_Ehrenfest_phase_average

