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
  complex(8) :: zCt_PT(2,2),zACt_PT(2)
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)
  complex(8) :: zSm(2,2),zHm(2,2),zHk(2,2),zHph(2,2),zHcoup(2,2),zDm(2,2)

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0; norm_t_l =0d0

  do itraj = 1,Ntraj
    if(myrank == 0)write(*,*)"itraj=",itraj
    call PTEF_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
    zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
    zphase = exp(-zI * aimag(sum(z_HO*conjg(zp_HO))))*(4d0/3d0)**Lsite


    call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
    zCt_PT(:,:) = 0d0; zCt_PT(1,1)=1d0; zCt_PT(2,2)=1d0

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      Xp_HO_new = 2d0*Xp_HO - Xp_HO_old + Fp_HO/mass*dt**2
      Vp_HO = 0.5d0*(Xp_HO_new - Xp_HO_old)/dt
!zs
      zACt_PT(:) = matmul(zSm, zCt_PT(:,2))
      zs = sum(conjg(zCt_PT(:,1))*zACt_PT(:))
      znorm_t = zs !1d0
!Ekin
      zACt_PT(:) = matmul(zHk, zCt_PT(:,2))
      zEkin_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs
!Eph
      zACt_PT(:) = matmul(zHph, zCt_PT(:,2))
      zEph_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs
!Eph
      zACt_PT(:) = matmul(zHcoup, zCt_PT(:,2))
      zEcoup_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs

      Ekin_l(it)=Ekin_l(it)+zEkin_t*zphase
      Eph_l(it)=Eph_l(it)+zEph_t*zphase
      Ecoup_l(it)=Ecoup_l(it)+zEcoup_t*zphase
      norm_t_l(it)=norm_t_l(it)+znorm_t*zphase


      call propagate_2x2_wavefunction(zCt_PT,zSm,zHm,zDm,dt*0.5d0)
      call dt_evolve_elec_pair
      X_HO_old = X_HO; X_HO = X_HO_new
      Xp_HO_old = Xp_HO; Xp_HO = Xp_HO_new
      call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
      call propagate_2x2_wavefunction(zCt_PT,zSm,zHm,zDm,dt*0.5d0)
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


subroutine propagate_2x2_wavefunction(zpsi,zSm,zHm,zDm,dt)
  implicit none
  complex(8),parameter :: zI = (0d0, 1d0)
  complex(8),intent(inout) :: zpsi(2,2)
  complex(8),intent(in) :: zSm(2,2),zHm(2,2),zDm(2,2)
  real(8),intent(in) :: dt
  complex(8) :: zSm_inv(2,2), zHeff(2,2)
  complex(8) :: ztpsi(2,2),zhpsi(2,2)
  real(8) :: det
  integer,parameter :: Ntaylor = 6
  integer :: itaylor
  complex(8) :: zfact

  det = zSm(1,1)*zSm(2,2)-zSm(1,2)*zSm(2,1)
  zSm_inv(1,1) = zSm(2,2)
  zSm_inv(2,2) = zSm(1,1)
  zSm_inv(1,2) = - zSm(1,2)
  zSm_inv(2,1) = - zSm(2,1)
  zSm_inv = zSm_inv/det

  zHeff = matmul(zSm_inv,(zHm - zDm))

  zfact = 1d0
  ztpsi = zpsi
  do itaylor = 1,Ntaylor
    zfact = zfact*(-zI*dt)/itaylor
    zhpsi = matmul(zHeff,ztpsi)
    zpsi = zpsi + zfact*zhpsi
    ztpsi = zhpsi
  end do

end subroutine propagate_2x2_wavefunction
