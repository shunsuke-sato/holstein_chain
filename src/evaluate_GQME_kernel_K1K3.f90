!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine evaluate_GQME_kernel_K1K3
  use global_variables
  implicit none
  complex(8),allocatable :: zK1_tmp(:,:,:),zK1_tmp_l(:,:,:)
  complex(8),allocatable :: zK3_tmp(:,:,:),zK3_tmp_l(:,:,:)
  real(8) :: phi1,phi2,phi3,phi4
  complex(8) :: zfact
  integer :: itraj,it,jsite
  integer :: a1,a2
  real(8),parameter :: rate = 0.1d0
  integer,parameter :: Nphase = 2
  integer  :: iphase,jphase
  real(8) :: phi(Nphase)

  allocate(zK1_tmp(Lsite,Lsite,0:Nt),zK1_tmp_l(Lsite,Lsite,0:Nt))
  allocate(zK3_tmp(Lsite,Lsite,0:Nt),zK3_tmp_l(Lsite,Lsite,0:Nt))


  zK1_tmp_l = 0d0; zK3_tmp_l = 0d0
! diagonal (1,1)
  do itraj = 1,Ntraj
    if(mod(itraj,Nprocs) /= myrank)cycle

    zC = 0d0; zC(1) = 1d0
!    call set_initial_conditions_elec
    call set_initial_conditions_ph(itraj)
    call calc_force_HO
    X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

    zfact = +zI*beta_KB*V_HO(1)

! == Kernel calc (1,1) ==
    do a1=1,Lsite
    do a2=1,Lsite
      zK1_tmp_l(a1,a2,0) = zK1_tmp_l(a1,a2,0) &
        +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
        *(X_HO(a1) - X_HO(a2))*zfact
      zK3_tmp_l(a1,a2,0) = zK3_tmp_l(a1,a2,0) &
        -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
    end do
    end do
! == Kernel calc (1,1) ==

    do it = 0,Nt-1

      X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
      V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt + F_HO/mass*dt

      call dt_evolve_elec
      X_HO_old = X_HO; X_HO = X_HO_new
      call calc_force_HO
! == Kernel calc (1,1) ==
      do a1=1,Lsite
      do a2=1,Lsite
         zK1_tmp_l(a1,a2,it+1) = zK1_tmp_l(a1,a2,it+1) &
              +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
              *(X_HO(a1) - X_HO(a2))*zfact
         zK3_tmp_l(a1,a2,it+1) = zK3_tmp_l(a1,a2,it+1) &
              -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
      end do
      end do
! == Kernel calc (1,1) ==

    end do

  end do

  call MPI_ALLREDUCE(zK1_tmp_l,zK1_tmp,(Nt+1)*Lsite**2, &
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(zK3_tmp_l,zK3_tmp,(Nt+1)*Lsite**2 &
    ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  zK1_tmp = zK1_tmp/dble(Ntraj);   zK3_tmp = zK3_tmp/dble(Ntraj)
!  zK1_tmp = zK1_tmp/dble(Ntraj)*Lsite**2;   zK3_tmp = zK3_tmp/dble(Ntraj)*Lsite**2
  zK1(:,:,1,1,:) = zK1_tmp(:,:,:); zK3(:,:,1,1,:) = zK3_tmp(:,:,:)
  if(myrank==0)  write(*,*)sum(abs(zK1(1,2,1,1,:))),sum(abs(zK1(1,2,1,2,:)))

! off-diagonal (1,j)
  do jsite = 2, Lsite
    call set_thermal_ph_dist
    zK1_tmp_l = 0d0; zK3_tmp_l = 0d0
    do iphase = 1,Nphase
       call random_number(phi(iphase)); phi(iphase) = 2d0 * pi *phi(iphase)
    end do
    do itraj = 1,Ntraj
      if(mod(itraj,Nprocs) /= myrank)cycle

      do iphase = 1,Nphase
      do jphase = 1,Nphase
! positive summ
!      call set_initial_conditions_elec
      zC = 0d0; zC(1) = sqrt(0.5d0); zC(jsite) = exp(zI*(phi(iphase) + 2d0*pi*dble(jphase)/dble(Nphase)  ))*sqrt(0.5d0)
!      zC(1) = exp(zI*phi1)*zC(1)
!      zC(jsite) = exp(zI*phi2)*zC(jsite)
      call set_initial_conditions_ph(itraj)
      call calc_force_HO
      X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2

      zfact = (X_HO(1) - X_HO(jsite))  +zI*0.5d0*beta_KB*(V_HO(1)+V_HO(jsite))
      zfact = zfact*zC(1)*conjg(zC(jsite))
      
! == Kernel calc (1,jsite) ==
      do a1=1,Lsite
      do a2=1,Lsite
        zK1_tmp_l(a1,a2,0) = zK1_tmp_l(a1,a2,0) &
          +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
          *(X_HO(a1) - X_HO(a2))*zfact
        zK3_tmp_l(a1,a2,0) = zK3_tmp_l(a1,a2,0) &
          -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
       end do
       end do
! == Kernel calc (1,jsite) ==

       do it = 0,Nt-1

         X_HO_new = 2d0*X_HO - X_HO_old + F_HO/mass*dt**2
         V_HO = 0.5d0*(X_HO_new - X_HO_old)/dt + F_HO/mass*dt

         call dt_evolve_elec
         X_HO_old = X_HO; X_HO = X_HO_new
         call calc_force_HO
! == Kernel calc (1,jsite) ==
         do a1=1,Lsite
         do a2=1,Lsite
            zK1_tmp_l(a1,a2,it+1) = zK1_tmp_l(a1,a2,it+1) &
                 +gamma**2*2d0*mass*omega0*zC(a2)*conjg(zC(a1)) &
                 *(X_HO(a1) - X_HO(a2))*zfact
            zK3_tmp_l(a1,a2,it+1) = zK3_tmp_l(a1,a2,it+1) &
                 -gamma*sqrt(2d0*mass*omega0)*zC(a2)*conjg(zC(a1))*zfact
         end do
         end do
! == Kernel calc (1,jsite) ==
       end do

    end do

    end do
    end do

      call MPI_ALLREDUCE(zK1_tmp_l,zK1_tmp,(Nt+1)*Lsite**2, &
        MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zK3_tmp_l,zK3_tmp,(Nt+1)*Lsite**2 &
        ,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      zK1_tmp = zK1_tmp/dble(Nphase**2*Ntraj)*2**2;   zK3_tmp = zK3_tmp/dble(Nphase**2*Ntraj)*2**2
      zK1(:,:,1,jsite,:) = zK1_tmp(:,:,:); zK3(:,:,1,jsite,:) = zK3_tmp(:,:,:)
  end do


end subroutine evaluate_GQME_kernel_K1K3
