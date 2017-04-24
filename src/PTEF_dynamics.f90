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
  complex(8) :: zSm(2,2),zVm(2,2),zVm_i_v22(2,2)
  complex(8) :: zHm(2,2),zHm_v(2,2)
  complex(8) :: zHk_v(2,2),zHph_v(2,2),zHcoup_v(2,2),zSm_v(2,2),zDm(2,2)
  complex(8) :: zHm_v_old(2,2)

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
    zHm_v_old = zHm_v
    zCt_PT(:,:) = zVm

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt
      Xp_HO_new = 2d0*Xp_HO - Xp_HO_old + Fp_HO/mass*dt**2
      Vp_HO = 0.5d0*(Xp_HO_new - Xp_HO_old)/dt
!zs
      zACt_PT(:) = matmul(zSm_v, zCt_PT(:,2))
      zs = sum(conjg(zCt_PT(:,1))*zACt_PT(:))
      znorm_t = zs !1d0
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
      zHm_v_old = zHm_v
      call PTEF_calc_matrix_old(zSm,zVm,zVm_i_v22,zHm,zHm_v,zHk_v,zHph_v,zHcoup_v,zSm_v,zDm)
      call propagate_2x2_wavefunction(zCt_PT,zHm_v,zHm_v_old,dt)
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


subroutine propagate_2x2_wavefunction(zpsi,zHm_v,zHm_v_old,dt)
  implicit none
  complex(8),parameter :: zI = (0d0, 1d0)
  complex(8),intent(inout) :: zpsi(2,2)
  complex(8),intent(in) :: zHm_v(2,2),zHm_v_old(2,2)
  real(8),intent(in) :: dt
  real(8) :: Am(2,2),Am_inv(2,2),b(2),x(2)
  complex(8) :: zHm(2,2),zFm(2,2)
  real(8) :: det, ss
  real(8) :: lambda(2),deps,eps_av
  complex(8) :: zeigv(2,2)
  complex(8) :: alpha,beta,zx,zy

  zHm = 0.5d0*(zHm_v + zHm_v_old)
  zFm = (zHm_v - zHm_v_old)/dt
  write(*,*)1,real(zFm(1,1)),real(zFm(2,2))

!  Am(1,1) = -2d0; Am(1,2) = (zHm(1,1)-zHm(2,2))*real(conjg(zHm(1,2))/zHm(1,2))
!  Am(2,1) = 0d0; Am(2,2) = (zHm(1,1)-zHm(2,2))*aimag(conjg(zHm(1,2))/zHm(1,2))
!  b(1) = real(zI*zFm(1,2)/zHm(1,2))
!  b(2) = aimag(zI*zFm(1,2)/zHm(1,2))

!  det = Am(1,1)*Am(2,2)-Am(1,2)*Am(2,1)
!  Am_inv(1,1) = Am(2,2)
!  Am_inv(2,2) = Am(1,1)
!  Am_inv(1,2) = - Am(1,2)
!  Am_inv(2,1) = - Am(2,1)
!  Am_inv = Am_inv/det

!  write(*,*)1,sum(Am_inv(1,:)*Am(:,1)),sum(Am_inv(1,:)*Am(:,2))
!  write(*,*)2,sum(Am_inv(2,:)*Am(:,1)),sum(Am_inv(2,:)*Am(:,2))

!  x = matmul(Am_inv,b)

!  alpha = x(1); beta =x(2)*conjg(zHm(1,2))
  alpha = real(zI*zFm(1,2)/zHm(1,2))/(zHm(1,1) - zHm(2,2))
  zHm(1,2) = zHm(1,2)*(1d0 + alpha); zHm(2,1) = conjg(zHm(1,2))
  

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
  

end subroutine propagate_2x2_wavefunction
