!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PTEF_dynamics
  use global_variables
  implicit none
  integer :: itraj,it
  complex(8) :: zEkin_t,zEph_t,zEcoup_t,znorm_t
  complex(8) :: zphase, zexponent0, zexponent,zs
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0; norm_t_l =0d0

  do itraj = 1,Ntraj
    if(myrank == 0)write(*,*)"itraj=",itraj
    call PTEF_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
    zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
    zphase = exp(-zI * aimag(sum(z_HO*conjg(zp_HO))))
    zexponent0 = -0.5d0*sum(abs(z_HO)**2 + abs(zp_HO)**2) +sum(conjg(z_HO)*zp_HO)

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      Xp_HO_new = 2d0*Xp_HO - Xp_HO_old + Fp_HO/mass*dt**2
      Vp_HO = 0.5d0*(Xp_HO_new - Xp_HO_old)/dt

      call calc_energy_pair(zEkin_t,zEph_t,zEcoup_t,znorm_t)
      z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
      zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
      zexponent = -0.5d0*sum(abs(z_HO)**2 + abs(zp_HO)**2) +sum(conjg(z_HO)*zp_HO)
!      zs = zphase*exp(zexponent - zexponent0) * (4d0/3d0)**Lsite
      zs = zphase * (4d0/3d0)**Lsite

      Ekin_l(it)=Ekin_l(it)+zEkin_t*zs
      Eph_l(it)=Eph_l(it)+zEph_t*zs
      Ecoup_l(it)=Ecoup_l(it)+zEcoup_t*zs
      norm_t_l(it)=norm_t_l(it)+znorm_t*zs

      call dt_evolve_elec_pair

      X_HO_old = X_HO; X_HO = X_HO_new
      Xp_HO_old = Xp_HO; Xp_HO = Xp_HO_new
      call calc_force_HO_pair

    end do

  end do

  call MPI_ALLREDUCE(Ekin_l,Ekin,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Eph_l,Eph,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(Ecoup_l,Ecoup,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(norm_t_l,norm_t,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  Ekin=Ekin/dble(Ntraj)
  Eph=Eph/dble(Ntraj)
  Ecoup=Ecoup/dble(Ntraj)
  norm_t=norm_t/dble(Ntraj)

  if(myrank == 0)then
    open(21,file="energy_t.out")
    write(21,"(A)")"# tt, Etot, Ekin, Eph, Ecoup"
    do it = 0,Nt
      write(21,"(999e26.16e3)")dt*it,Ekin(it)+Eph(it)+Ecoup(it)&
           ,Ekin(it),Eph(it),Ecoup(it),norm_t(it)
    end do
    close(21)
  end if

end subroutine PTEF_dynamics
