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
!  complex(8) :: zSm(2,2),zHm(2,2),zHk(2,2),zHph(2,2),zHcoup(2,2),zDm(2,2)
  complex(8) :: zSm(2,2),zVm(2,2),zHm(2,2),zDm(2,2)
  complex(8) :: zSm_old(2,2),zDm_old(2,2),zHm_old(2,2)
  complex(8) :: zHk(2,2),zHph(2,2),zHcoup(2,2)
  complex(8) :: znorm0
  complex(8) :: zdrho(Lsite,Lsite)
  integer :: isite,jsite,ialpha

  do isite = 1,Lsite
    do jsite = 1,Lsite
      zdrho(isite,jsite) = &
        exp(zI*2d0*pi*(isite-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite)) &
        *exp(-zI*2d0*pi*(jsite-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite))
    end do
  end do

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0; norm_t_l =0d0

  do itraj = 1,Ntraj
    if(myrank == 0)write(*,*)"itraj=",itraj
    ialpha = mod(itraj,Lsite) + 1
    call PTEF_initial_condition(ialpha)
    if(mod(itraj,Nprocs) /= myrank)cycle
    z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
    zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
!    zphase = exp(-zI * aimag(sum(z_HO*conjg(zp_HO))))*(4d0/3d0)**Lsite
    zphase = exp(0.5d0*sum(abs(z_HO - zp_HO)**2))*(4d0/3d0)**Lsite
    zphase = zphase * zdrho(1,ialpha)*Lsite**2


    call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
    zCt_PT(:,:) = 0d0; zCt_PT(1,1)=1d0; zCt_PT(2,2)=1d0
    zACt_PT(:) = matmul(zSm, zCt_PT(:,2))
!    znorm0 = sum(conjg(zCt_PT(:,1))*zACt_PT(:))


    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt + F_HO/mass*dt
!      V_HO = V_HO+ F_HO/mass*dt
      Xp_HO_new = 2d0*Xp_HO - Xp_HO_old + Fp_HO/mass*dt**2
      Vp_HO = 0.5d0*(Xp_HO_new - Xp_HO_old)/dt + Fp_HO/mass*dt
!      Vp_HO = Vp_HO + Fp_HO/mass*dt
!zs
      zACt_PT(:) = matmul(zSm, zCt_PT(:,2))
      zs = sum(conjg(zCt_PT(:,1))*zACt_PT(:))
!      znorm_t = zs/znorm0
      znorm_t = zs
!      znorm_t = (zHm(1,1)+zHm(2,2))/znorm0
!      znorm_t = (zHm(1,1)+zHm(2,2))/znorm0

!Ekin
      zACt_PT(:) = matmul(zHk, zCt_PT(:,2))
!      zEkin_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/znorm0
      zEkin_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))
!Eph
      zACt_PT(:) = matmul(zHph, zCt_PT(:,2))
!      zEph_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/znorm0
      zEph_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))
!Eph
      zACt_PT(:) = matmul(zHcoup, zCt_PT(:,2))
!      zEcoup_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/znorm0
      zEcoup_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))

      Ekin_l(it)=Ekin_l(it)+zEkin_t*zphase
      Eph_l(it)=Eph_l(it)+zEph_t*zphase
      Ecoup_l(it)=Ecoup_l(it)+zEcoup_t*zphase
      norm_t_l(it)=norm_t_l(it)+znorm_t*zphase


      call dt_evolve_elec_pair
      X_HO_old = X_HO; X_HO = X_HO_new
      Xp_HO_old = Xp_HO; Xp_HO = Xp_HO_new
!      call PTET_dt_evolve(zCt_PT,zSm,zHm,zDm,dt*0.5d0)
      zSm_old=zSm; zDm_old=zDm; zHm_old = zHm
      call calc_force_HO_pair
      call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
      call PTET_dt_evolve(zCt_PT,zSm,zHm,zDm,zSm_old,zHm_old,zDm_old,dt)
!      call PTEF_calc_matrix_old(zSm,zVm,zVm_i_v22,zHm,zHm_v,zHk_v,zHph_v,zHcoup_v,zSm_v,zDm)
!      call propagate_2x2_wavefunction(zCt_PT,zHm_v,zHm_v_old,zSm,zSm_old,dt)


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
subroutine PTET_dt_evolve(zpsi,zSm,zHm,zDm,zSm_old,zHm_old,zDm_old,dt)
  implicit none
  complex(8),parameter :: zI = (0d0,1d0)
  complex(8),intent(inout) :: zpsi(2,2)
  complex(8),intent(in) :: zSm(2,2),zHm(2,2),zDm(2,2)
  complex(8),intent(in) :: zSm_old(2,2),zHm_old(2,2),zDm_old(2,2)
  real(8),intent(in) :: dt
  complex(8) :: zHeff(2,2),zHeff_old(2,2),zProp(2,2)
  complex(8) :: zAm(2,2), zAm_inv(2,2), zfact
  complex(8) :: ztpsi(2,2),zhtpsi(2,2)
  integer :: i
  integer,parameter :: NTex = 6
  integer :: nPropagator
  integer,parameter :: num_Taylor = 0
  integer,parameter :: num_CN = 1


  nPropagator = num_Taylor
!  nPropagator = num_CN

  select case(nPropagator)

  case(num_CN)
    zAm = zSm
    call zinv2x2mat(zAm,zAm_inv) 
    zHeff = matmul(zAm_inv,(zHm-zDm))

    zAm = zSm_old
    call zinv2x2mat(zAm,zAm_inv) 
    zHeff_old = matmul(zAm_inv,(zHm_old-zDm_old))
    
    zAm = 0d0; zAm(1,1) = 1d0; zAm(2,2) = 1d0
    zAm = zAm + zI*0.5d0*dt*zHeff
    call zinv2x2mat(zAm,zAm_inv) 
    
    zAm = 0d0; zAm(1,1) = 1d0; zAm(2,2) = 1d0
    zAm = zAm - zI*0.5d0*dt*zHeff_old
    zProp = matmul(zAm_inv,zAm)
    
    zpsi = matmul(zProp,zpsi)

  case(num_Taylor)

    zAm = zSm_old
    call zinv2x2mat(zAm,zAm_inv) 
    zHeff_old = matmul(zAm_inv,(zHm_old-zDm_old))

    zfact = 1d0
    ztpsi = zpsi
    do i = 1,NTex
      zfact = zfact*(-zI*dt*0.5d0)/i
      zhtpsi = matmul(zHeff_old,ztpsi)
      zpsi = zpsi + zfact*zhtpsi
      ztpsi = zhtpsi
    end do

    zAm = zSm
    call zinv2x2mat(zAm,zAm_inv) 
    zHeff = matmul(zAm_inv,(zHm-zDm))

    zfact = 1d0
    ztpsi = zpsi
    do i = 1,NTex
      zfact = zfact*(-zI*dt*0.5d0)/i
      zhtpsi = matmul(zHeff,ztpsi)
      zpsi = zpsi + zfact*zhtpsi
      ztpsi = zhtpsi
    end do

  case default 
    stop 'Invalid propagator'
  end select

end subroutine PTET_dt_evolve


subroutine PTEF_dynamics_old
  use global_variables
  implicit none
  integer :: itraj,it
  complex(8) :: zEkin_t,zEph_t,zEcoup_t,znorm_t
  complex(8) :: zphase, zexponent0, zexponent,zs
  complex(8) :: zCt_PT(2,2),zACt_PT(2)
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)
!  complex(8) :: zSm(2,2),zHm(2,2),zHk(2,2),zHph(2,2),zHcoup(2,2),zDm(2,2)
  complex(8) :: zSm(2,2),zVm(2,2),zVm_i_v22(2,2)
  complex(8) :: zHm(2,2),zHm_v(2,2)
  complex(8) :: zHk_v(2,2),zHph_v(2,2),zHcoup_v(2,2),zSm_v(2,2),zDm(2,2)
  complex(8) :: zHm_v_old(2,2),zSm_old(2,2)

  Ekin_l = 0d0; Eph_l = 0d0; Ecoup_l = 0d0; norm_t_l =0d0

  do itraj = 1,Ntraj
    if(myrank == 0)write(*,*)"itraj=",itraj
    call PTEF_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
    zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
    zphase = exp(-zI * aimag(sum(z_HO*conjg(zp_HO))))*(4d0/3d0)**Lsite


!    call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
!    zCt_PT(:,:) = 0d0; zCt_PT(1,1)=1d0; zCt_PT(2,2)=1d0
    call PTEF_calc_matrix_old(zSm,zVm,zVm_i_v22,zHm,zHm_v,zHk_v,zHph_v,zHcoup_v,zSm_v,zDm)
    zHm_v_old = zHm_v; zSm_old = zSm
    zCt_PT(:,:) = zVm

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      Xp_HO_new = 2d0*Xp_HO - Xp_HO_old + Fp_HO/mass*dt**2
      Vp_HO = 0.5d0*(Xp_HO_new - Xp_HO_old)/dt
!zs
      zACt_PT(:) = matmul(zSm_v, zCt_PT(:,2))
      zs = sum(conjg(zCt_PT(:,1))*zACt_PT(:))

      zACt_PT(:) = matmul(zHm_v, zCt_PT(:,2))
      znorm_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:)) !/zs
!      znorm_t = zs !1d0
!Ekin
      zACt_PT(:) = matmul(zHk_v, zCt_PT(:,2))
      zEkin_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs
!Eph
      zACt_PT(:) = matmul(zHph_v, zCt_PT(:,2))
      zEph_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs
!Eph
      zACt_PT(:) = matmul(zHcoup_v, zCt_PT(:,2))
      zEcoup_t = sum(conjg(zCt_PT(:,1))*zACt_PT(:))/zs

      Ekin_l(it)=Ekin_l(it)+zEkin_t*zphase
      Eph_l(it)=Eph_l(it)+zEph_t*zphase
      Ecoup_l(it)=Ecoup_l(it)+zEcoup_t*zphase
      norm_t_l(it)=norm_t_l(it)+znorm_t*zphase


      call dt_evolve_elec_pair
      X_HO_old = X_HO; X_HO = X_HO_new
      Xp_HO_old = Xp_HO; Xp_HO = Xp_HO_new
!      call PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
      zHm_v_old = zHm_v; zSm_old = zSm
      call PTEF_calc_matrix_old(zSm,zVm,zVm_i_v22,zHm,zHm_v,zHk_v,zHph_v,zHcoup_v,zSm_v,zDm)
      call propagate_2x2_wavefunction(zCt_PT,zHm_v,zHm_v_old,zSm,zSm_old,dt)
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

end subroutine PTEF_dynamics_old


subroutine propagate_2x2_wavefunction(zpsi,zHm_v,zHm_v_old,zSm,zSm_old,dt)
  implicit none
  complex(8),parameter :: zI = (0d0, 1d0)
  complex(8),intent(inout) :: zpsi(2,2)
  complex(8),intent(in) :: zHm_v(2,2),zHm_v_old(2,2)
  complex(8),intent(in) :: zSm(2,2),zSm_old(2,2)
  real(8),intent(in) :: dt
  real(8) :: Am(2,2),Am_inv(2,2),b(2),x(2)
  complex(8) :: zHm_eff(2,2),zFm_eff(2,2),zSm_eff(2,2)
  complex(8) :: zVint_eff(2,2)
  complex(8) :: zVm(2,2),zVm_i_v22(2,2)
  complex(8) :: zVm_old(2,2),zVm_i_v22_old(2,2)
  complex(8) :: zF_val
  real(8) :: det, ss, alpha
  complex(8) :: zpsi_old(2,2),zvec(2), zMm1(2,2), zMm2(2,2), zs
  complex(8) :: zAm(2,2),zAm_inv(2,2)
  complex(8) :: zAop(2,2),zBop(2,2),zCop(2,2)
  
!  return
  zHm_eff = 0.5d0*(zHm_v + zHm_v_old)
  zFm_eff = -zI*(zHm_v - zHm_v_old)/dt
  zpsi_old = zpsi

  zAm = zHm_eff
  call zinv2x2mat(zAm,zAm_inv) 
  zMm1 = (zHm_v - zHm_v_old)/dt
  zMm2 = matmul(zMm1,zAm_inv);  zAop = -zI*0.5d0* ( matmul(zAm_inv,zMm1) - matmul(zMm1,zAm_inv) )
!  zBop = zMm1
!  zMm2 = matmul(zMm1,zAm_inv);  zBop = matmul(zMm2,zMm1)
!  zMm2 = matmul(zAm_inv,zMm1);  zBop = 0.5d0*( zBop + matmul(zAm_inv,zMm2))

!  write(*,*)"a11",zAm(1,1)
!  write(*,*)"a12",zAm(1,2)
!  write(*,*)"a21",zAm(2,1)
!  write(*,*)"a22",zAm(2,2)

! == First step; Predictor for the effective interction zVint
  zvec = matmul(zFm_eff,zpsi(:,2))
  zF_val = sum(conjg(zpsi(:,1))*zvec)

! first; Aop
!  zMm1 = -zi*( matmul(zHm_v_old,zAop) - matmul(zAop,zHm_v_old) )
  zMm2 = matmul(zHm_v_old,zAop) - matmul(zAop,zHm_v_old) 
  zvec = matmul(zMm2,zpsi(:,2))
  zs = sum(conjg(zpsi(:,1))*zvec)
  alpha = real(zF_val/zs)


!  write(*,*)'error',abs(alpha*zs - zF_val)
!  Am(1,1) = real(zs); Am(2,1) = aimag(zs)
!
!! second; Bop
!  zMm1 = -zi*( matmul(zHm_v_old,zBop) - matmul(zBop,zHm_v_old) )
!  zMm2 = matmul(zHm_v_old,zMm1) - matmul(zMm1,zHm_v_old) 
!  zvec = matmul(zMm2,zpsi(:,2))
!  zs = sum(conjg(zpsi(:,1))*zvec)
!  Am(1,2) = real(zs); Am(2,2) = aimag(zs)
!
!  call inv2x2mat(Am,Am_inv)
!!  write(*,*)sum(Am(1,:)*Am_inv(:,1)),sum(Am(1,:)*Am_inv(:,2))
!  b(1) = real(zF_val); b(2) = aimag(zF_val)
!  x = matmul(Am_inv,b)
!
!  zCop = x(1)*zAop + x(2)*zBop
!!  write(*,*)x
  alpha = 1d0
  zHm_eff = zHm_v_old + alpha * zAop 
  call prop(zpsi,zHm_eff,dt*0.5d0)
  zHm_eff = zHm_v + alpha * zAop  
  call prop(zpsi,zHm_eff,dt*0.5d0)

  return ! return for check
! == Second step; Predictor for the effective interction zVint


  zVm(1,1) = 1d0; zVm(1,2) = zSm(1,2)
  zVm(2,1) = 0d0; zVm(2,2) = sqrt(abs(1d0 - abs(zSm(1,2))**2))
  zVm_i_v22(1,1) = zVm(2,2); zVm_i_v22(1,2) = - zVm(1,2)
  zVm_i_v22(2,1) = 0d0     ; zVm_i_v22(2,2) =   zVm(1,1)

  zvec = matmul(zFm_eff,zpsi(:,2))
  zF_val = 0.5d0*( zF_val + sum(conjg(zpsi(:,1))*zvec) )

! first
  zMm1(1,1) = 0d0; zMm1(1,2) = 1d0
  zMm1(2,1) = 1d0; zMm1(2,2) = 0d0
  zMm2 = matmul(zMm1,zVm_i_v22); zMm1 = matmul( transpose(conjg(zVm_i_v22)),zMm2 )
  zMm1 = zMm1/abs(zVm(2,2))**2
  zMm2 = matmul(zHm_v,zMm1) - matmul(zMm1,zHm_v)
  zvec = matmul(zMm2,zpsi(:,2))
  zs = sum(conjg(zpsi(:,1))*zvec)
  Am(1,1) = Am(1,1) + real(zs); Am(2,1) = Am(2,1) + aimag(zs)

! second
  zMm1(1,1) = 0d0; zMm1(1,2) = zI
  zMm1(2,1) = -zI; zMm1(2,2) = 0d0
  zMm2 = matmul(zMm1,zVm_i_v22); zMm1 = matmul( transpose(conjg(zVm_i_v22)),zMm2 )
  zMm1 = zMm1/abs(zVm(2,2))**2
  zMm2 = matmul(zHm_v,zMm1) - matmul(zMm1,zHm_v)
  zvec = matmul(zMm2,zpsi(:,2))
  zs = sum(conjg(zpsi(:,1))*zvec)
  Am(1,2) = Am(1,2) + real(zs); Am(2,2) = Am(2,2) + aimag(zs)
  Am = Am/2d0

  call inv2x2mat(Am,Am_inv)
  b(1) = real(zF_val); b(2) = aimag(zF_val)
  x = matmul(Am_inv,b)


! == Final step; True propagation

  zpsi = zpsi_old

  zMm1(1,1) = 0d0; zMm1(1,2) = x(1) + zI*x(2)
  zMm1(2,1) = x(1)-zI*x(2); zMm1(2,2) = 0d0
  zMm2 = matmul(zMm1,zVm_i_v22_old); zMm1 = matmul( transpose(conjg(zVm_i_v22_old)),zMm2 )
  zVint_eff = zMm1/abs(zVm_old(2,2))**2

  zHm_eff = zHm_v_old + zVint_eff
  call prop(zpsi,zHm_eff,dt*0.5d0)

  zMm1(1,1) = 0d0; zMm1(1,2) = x(1) + zI*x(2)
  zMm1(2,1) = x(1)-zI*x(2); zMm1(2,2) = 0d0
  zMm2 = matmul(zMm1,zVm_i_v22); zMm1 = matmul( transpose(conjg(zVm_i_v22)),zMm2 )
  zVint_eff = zMm1/abs(zVm(2,2))**2

  zHm_eff = zHm_v + zVint_eff
  call prop(zpsi,zHm_eff,dt*0.5d0)

end subroutine propagate_2x2_wavefunction

subroutine inv2x2mat(Am,Am_inv)
  implicit none
  real(8),intent(in) :: Am(2,2)
  real(8),intent(out) :: Am_inv(2,2)
  real(8) :: det

  det = Am(1,1)*Am(2,2)-Am(1,2)*Am(2,1)
!  write(*,*)'det',det
  Am_inv(1,1) = Am(2,2)
  Am_inv(2,2) = Am(1,1)
  Am_inv(1,2) = - Am(1,2)
  Am_inv(2,1) = - Am(2,1)
  Am_inv = Am_inv/det

end subroutine inv2x2mat

subroutine zinv2x2mat(zAm,zAm_inv)
  implicit none
  complex(8),intent(in) :: zAm(2,2)
  complex(8),intent(out) :: zAm_inv(2,2)
  complex(8) :: zdet

  zdet = zAm(1,1)*zAm(2,2)-zAm(1,2)*zAm(2,1)
!  write(*,*)'zdet',zdet
  zAm_inv(1,1) = zAm(2,2)
  zAm_inv(2,2) = zAm(1,1)
  zAm_inv(1,2) = - zAm(1,2)
  zAm_inv(2,1) = - zAm(2,1)
  zAm_inv = zAm_inv/zdet

end subroutine zinv2x2mat

subroutine prop(zpsi,zHm,dt)
  implicit none
  complex(8),parameter :: zI = (0d0, 1d0)
  complex(8),intent(inout) :: zpsi(2,2)
  complex(8),intent(in) :: zHm(2,2)
  real(8),intent(in) :: dt
  real(8) :: lambda(2),deps,eps_av, ss
  complex(8) :: zeigv(2,2)
  complex(8) :: alpha,beta,zx,zy


  eps_av = 0.5d0*( zHm(1,1) + zHm(2,2) )
  deps = sqrt( abs(zHm(1,1) - zHm(2,2))**2 + 4d0*abs(zHm(1,2))**2   )
  if(real(zHm(1,1) - zHm(2,2) ) > 0d0 )then
     lambda(1) = eps_av + 0.5d0*deps
     lambda(2) = eps_av - 0.5d0*deps
  else
     lambda(1) = eps_av - 0.5d0*deps
     lambda(2) = eps_av + 0.5d0*deps
  end if

  zy = conjg(zHm(1,2))/(lambda(1) - zHm(2,2))
  ss = 1d0 + abs(zy)**2; ss = 1d0/sqrt(ss)
  zeigv(1,1) = ss; zeigv(2,1) = zy*ss

  zx = zHm(1,2)/(lambda(2) - zHm(1,1))
  ss = 1d0 + abs(zx)**2; ss = 1d0/sqrt(ss)
  zeigv(1,2) = zx*ss; zeigv(2,2) = ss

! psi 1
  zx = sum(conjg(zeigv(:,1))*zpsi(:,1))*exp(-zI*dt*lambda(1))
  zy = sum(conjg(zeigv(:,2))*zpsi(:,1))*exp(-zI*dt*lambda(2))
  zpsi(:,1) = zx * zeigv(:,1) + zy * zeigv(:,2)

! psi 2
  zx = sum(conjg(zeigv(:,1))*zpsi(:,2))*exp(-zI*dt*lambda(1))
  zy = sum(conjg(zeigv(:,2))*zpsi(:,2))*exp(-zI*dt*lambda(2))
  zpsi(:,2) = zx * zeigv(:,1) + zy * zeigv(:,2)
end subroutine prop
