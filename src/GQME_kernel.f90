!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GQME_kernel
  use global_variables
  implicit none
  complex(8) :: zK1_tmp(Lsite,Lsite),zK1_tmp_l(Lsite,Lsite)
  complex(8) :: zK3_tmp(Lsite,Lsite),zK3_tmp_l(Lsite,Lsite)
  real(8) :: phi, zfact
  integer :: itraj,it,jsite

  call allocate_GQME_kernel

  zK1_tmp_l = 0d0; zK3_tmp_l = 0d0
! diagonal (1,1)
  do itraj = 1,Ntraj
    if(mod(itraj,Nprocs) /= myrank)cycle

    zC = 0d0; zC(1) = 1d0
    call set_initial_conditions_ph(itraj)
    call calc_force_HO
    X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

    zfact = -zI*beta_KB*V_HO(1)

    do it = 0,Nt

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

      call dt_evolve_elec
      X_HO_old = X_HO; X_HO = X_HO_new
      call calc_force_HO
! Kernel calc
      do a1=1,Lsite
      do a2=1,Lsite
         zK1_tmp_l(a1,a2) = zK1_tmp_l(a1,a2) &
              +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
              *(X_HO(a1) - X_HO(a2))*zfact
         zK3_tmp_l(a1,a2) = zK3_tmp_l(a1,a2) &
              -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
      end do
      end do
    end do

  end do

  call MPI_ALLREDUCE(zK1_tmp_l,zK1_tmp,Lsite**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(zK3_tmp_l,zK3_tmp,Lsite**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  zK1_tmp = zK1_tmp/dble(Ntraj);   zK3_tmp = zK3_tmp/dble(Ntraj)
  zK1(:,:,1,1) = zK1_tmp(:,:); zK3(:,:,1,1) = zK3_tmp(:,:)

  zK1_tmp_l = 0d0; zK3_tmp_l = 0d0
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

      zfact = (X_HO(1) - X_HO(jsite))  -zI*0.5d0*beta_KB*(V_HO(1)+V_HO(jsite))
      zfact = zfact*zC(1)*conjg(zC(jsite))
      
      do it = 0,Nt

         X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
         V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

         call dt_evolve_elec
         X_HO_old = X_HO; X_HO = X_HO_new
         call calc_force_HO
! Kernel calc
         do a1=1,Lsite
         do a2=1,Lsite
            zK1_tmp_l(a1,a2) = zK1_tmp_l(a1,a2) &
                 +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
                 *(X_HO(a1) - X_HO(a2))*zfact
            zK3_tmp_l(a1,a2) = zK3_tmp_l(a1,a2) &
                 -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
         end do
         end do
      end do

! negative summ
      zC = 0d0; zC(1) = sqrt(0.5d0); zC(jsite) = -exp(zI*phi)*sqrt(0.5d0)
      call set_initial_conditions_ph(itraj)
      call calc_force_HO
      X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

      zfact = (X_HO(1) - X_HO(jsite))  -zI*0.5d0*beta_KB*(V_HO(1)+V_HO(jsite))
      zfact = zfact*zC(1)*conjg(zC(jsite))
      
      do it = 0,Nt

         X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
         V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt

         call dt_evolve_elec
         X_HO_old = X_HO; X_HO = X_HO_new
         call calc_force_HO
! Kernel calc
         do a1=1,Lsite
         do a2=1,Lsite
            zK1_tmp_l(a1,a2) = zK1_tmp_l(a1,a2) &
                 +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
                 *(X_HO(a1) - X_HO(a2))*zfact
            zK3_tmp_l(a1,a2) = zK3_tmp_l(a1,a2) &
                 -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
         end do
         end do
      end do

      call MPI_ALLREDUCE(zK1_tmp_l,zK1_tmp,Lsite**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zK3_tmp_l,zK3_tmp,Lsite**2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      zK1_tmp = zK1_tmp/dble(Ntraj)*2d0;   zK3_tmp = zK3_tmp/dble(Ntraj)*2d0
      zK1(:,:,1,jsite) = zK1_tmp(:,:); zK3(:,:,1,jsite) = zK3_tmp(:,:)
    end do
  end do


end subroutine GQME_kernel
