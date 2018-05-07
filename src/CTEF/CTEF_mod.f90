!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module CTEF_mod
  use global_variables
  use math_module
  implicit none

  private

  integer,parameter :: Nphase = 2, Nscf_refine = 2, Nscf_pred_corr = 2
  real(8),parameter :: epsilon_norm = 10d0/1d2 ! 10%
  real(8),parameter :: sigma_correlated_gaussian = 1.0d0 !1d0
  real(8),parameter :: sigma_independent_gaussian = sqrt(0.5d0) !1d0
  integer,parameter :: nsize_store = 10000
  complex(8),allocatable :: zpsi_CTEF(:,:),zHO_CTEF(:,:)
  complex(8),allocatable :: zHO_dot_CTEF(:,:)

  complex(8) :: zSs_CTEF(2,2), zDs_CTEF(2,2)
  complex(8) :: zSb_inv_CTEF(2,2)
  complex(8) :: zSb_CTEF(2,2), zDb_CTEF(2,2)
  complex(8) :: zEs_CTEF(2,2), zEc_CTEF(2,2), zEb_CTEF(2,2)
  complex(8) :: zHb_eff_CTEF(2,2), zSsb_CTEF(2,2), zSsb_inv_CTEF(2,2)
  complex(8),allocatable :: zX_HO_CTEF(:,:,:), zF_HO_CTEF(:,:)


  integer,parameter :: bath_propagator_direct = 0
  integer,parameter :: bath_propagator_taylor = 1
  integer,parameter :: bath_propagator_taylor_mod = 2
  integer,parameter :: bath_propagator_diag = 3
  integer,parameter :: bath_propagator_CrankNicolson = 4
  integer,parameter :: iflag_bath_propagator = bath_propagator_direct

  public :: CTEF, CTEF_kernel

  contains

!-----------------------------------------------------------------------------------------
    subroutine CTEF
      implicit none

      call CTEF_allocation
      call CTEF_dynamics

    end subroutine CTEF
!-----------------------------------------------------------------------------------------
    subroutine CTEF_kernel
      implicit none

      call CTEF_allocation
      call CTEF_dynamics_kernel

    end subroutine CTEF_kernel
!-----------------------------------------------------------------------------------------
    subroutine CTEF_allocation
      implicit none
      integer :: i,j

      allocate(zpsi_CTEF(Lsite,2),zHO_CTEF(Lsite,2),Hmat_kin(Lsite,Lsite))
      allocate(zHO_dot_CTEF(Lsite,2))
      allocate(zX_HO_CTEF(Lsite,2,2), zF_HO_CTEF(Lsite,2))

      Hmat_kin = 0d0
      do i =1,Lsite
        j=i+1
        if(j>Lsite)j=j-Lsite
        Hmat_kin(i,j) = -t0
        j=i-1
        if(j<1)j=j+Lsite
        Hmat_kin(i,j) = -t0
      end do

    end subroutine CTEF_allocation
!-----------------------------------------------------------------------------------------
    subroutine CTEF_dynamics
      implicit none
      complex(8) :: zpsi_store(Lsite,2),zHO_store(Lsite,2),zHO_gauss_store(Lsite,2)
      complex(8) :: zweight, zweight0
      real(8) :: norm, phi0,phi
      integer :: iphase, it, i,j
      integer :: itraj_t, istore
      integer(8) :: itraj
      integer :: i_antithetic, j_antithetic
      integer :: i_dm, j_dm
      real(8) :: norm_CTEF(0:Nt+1),norm_CTEF_l(0:Nt+1),norm_CTEF_t(0:Nt+1)
      real(8) :: norm_CTEF_sl(0:Nt+1)
      real(8) :: norm_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ekin_CTEF(0:Nt+1),Ekin_CTEF_l(0:Nt+1),Ekin_CTEF_t(0:Nt+1)
      real(8) :: Ekin_CTEF_sl(0:Nt+1)
      real(8) :: Ekin_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ebath_CTEF(0:Nt+1),Ebath_CTEF_l(0:Nt+1),Ebath_CTEF_t(0:Nt+1)
      real(8) :: Ebath_CTEF_sl(0:Nt+1)
      real(8) :: Ebath_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ecoup_CTEF(0:Nt+1),Ecoup_CTEF_l(0:Nt+1),Ecoup_CTEF_t(0:Nt+1)
      real(8) :: Ecoup_CTEF_sl(0:Nt+1)
      real(8) :: Ecoup_CTEF_phase_ave(0:Nt+1)
      complex(8) :: zrho_dm
      real(8) :: x1,x2,p1,p2
      integer,parameter :: ran_len = 1
      real(8) :: rvec(ran_len)
      logical :: is_norm_converged
      integer(8) :: ntraj_tot, ntraj_tot_l
      integer(8) :: ntraj_stable, ntraj_stable_l

      norm_CTEF_l = 0d0; Ekin_CTEF_l = 0d0
      Ebath_CTEF_l = 0d0; Ecoup_CTEF_l = 0d0
      ntraj_tot_l = 0; ntraj_stable_l = 0
      
      itraj = 0
      do itraj_t  = 1, max(Ntraj/nsize_store,1)
      norm_CTEF_sl = 0d0; Ekin_CTEF_sl = 0d0
      Ebath_CTEF_sl = 0d0; Ecoup_CTEF_sl = 0d0
      do istore = 1, min(nsize_store,Ntraj)
        itraj = itraj + 1
! == bath distribution
        call bath_sampling_correlated_gaussian(zHO_store,zweight0)
!        call bath_sampling_independent_gaussian(zHO_store,zweight0)

!        do i = 1,Lsite
!          call gaussian_random_number(x1,p1)
!          call gaussian_random_number(x2,p2)
!          zHO_gauss_store(i,1) = (x1 + zI * p1)*sqrt(2d0/3d0)
!          zHO_gauss_store(i,2) = (x2 + zI * p2)*sqrt(0.5d0)
!        end do
!        zHO_store(:,1) = zHO_gauss_store(:,1)
!        zHO_store(:,2) = zHO_gauss_store(:,2)+0.5d0*zHO_store(:,1)
!        call calc_zweight(zHO_store,zweight0)

! == bath distribution

        call ranlux_double (rvec, ran_len)
        phi0 = rvec(1); phi0 = 2d0*pi*phi0
        if(myrank == 0 .and. mod(itraj,max(Ntraj/200,1))==0)write(*,*)"itraj=",itraj,"/",Ntraj
        if(mod(itraj,Nprocs) /= myrank)cycle


! == localized init wf
        i_dm = 1
        do j_dm = 1,Lsite
          zpsi_store = 0d0 
          zpsi_store(i_dm,1) = 1d0; zpsi_store(j_dm,2) = 1d0 
          zrho_dm = exp(zI*2d0*pi*(j_dm-i_dm)*dble(Lsite/2)/dble(Lsite)) !&
              !              /dble(Lsite)*dble(Lsite*)
! == localized init wf

          zweight = zweight0 * zrho_dm

          norm_CTEF_phase_ave = 0d0
          Ekin_CTEF_phase_ave = 0d0
          Ebath_CTEF_phase_ave = 0d0
          Ecoup_CTEF_phase_ave = 0d0
          is_norm_converged = .true.
          do iphase = 1,Nphase
            if(.not. is_norm_converged)exit
            phi = phi0 + 2d0*pi*dble(iphase-1)/Nphase
            call set_initial_condition(zpsi_store,zHO_store, &
              zpsi_CTEF, zHO_CTEF, phi, norm)
            
            call propagation(norm_CTEF_t,Ekin_CTEF_t, Ebath_CTEF_t, Ecoup_CTEF_t)
              !          if(myrank == 0)write(*,*)"norm",norm_CTEF_t(0),norm
            norm_CTEF_phase_ave = norm_CTEF_phase_ave &
              + norm_CTEF_t*exp(-zI*phi)*norm*zweight
            Ekin_CTEF_phase_ave = Ekin_CTEF_phase_ave &
              + Ekin_CTEF_t*exp(-zI*phi)*norm*zweight
            Ebath_CTEF_phase_ave = Ebath_CTEF_phase_ave &
              + Ebath_CTEF_t*exp(-zI*phi)*norm*zweight
            Ecoup_CTEF_phase_ave = Ecoup_CTEF_phase_ave &
              + Ecoup_CTEF_t*exp(-zI*phi)*norm*zweight
            
            if(.not. abs(norm_CTEF_t(Nt)-1d0) < epsilon_norm )is_norm_converged = .false.
            
          end do
          if(is_norm_converged)then
            ntraj_stable_l = ntraj_stable_l + 1
            norm_CTEF_sl = norm_CTEF_sl + norm_CTEF_phase_ave
            Ekin_CTEF_sl = Ekin_CTEF_sl + Ekin_CTEF_phase_ave
            Ebath_CTEF_sl = Ebath_CTEF_sl + Ebath_CTEF_phase_ave
            Ecoup_CTEF_sl = Ecoup_CTEF_sl + Ecoup_CTEF_phase_ave
          end if
          ntraj_tot_l = ntraj_tot_l + 1

        end do

      end do

      norm_CTEF_l = norm_CTEF_l + norm_CTEF_sl/dble(nsize_store)
      Ekin_CTEF_l = Ekin_CTEF_l + Ekin_CTEF_sl/dble(nsize_store)
      Ebath_CTEF_l = Ebath_CTEF_l + Ebath_CTEF_sl/dble(nsize_store)
      Ecoup_CTEF_l = Ecoup_CTEF_l + Ecoup_CTEF_sl/dble(nsize_store)


      end do

      call MPI_ALLREDUCE(ntraj_tot_l,ntraj_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ntraj_stable_l,ntraj_stable,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)

      norm_CTEF_l = norm_CTEF_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)
      Ekin_CTEF_l = Ekin_CTEF_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)
      Ebath_CTEF_l = Ebath_CTEF_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)
      Ecoup_CTEF_l = Ecoup_CTEF_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)
      call MPI_ALLREDUCE(norm_CTEF_l,norm_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ekin_CTEF_l,Ekin_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ebath_CTEF_l,Ebath_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ecoup_CTEF_l,Ecoup_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myrank == 0)then
        write(*,*)"# of total   trajectries",ntraj_tot
        write(*,*)"# of stable  trajectries",ntraj_stable
        write(*,*)"# of skipped trajectries",ntraj_tot - ntraj_stable
        open(21,file="CTEF_norm.out")
        do it = 0,Nt
          write(21,"(999e26.16e3)")dt*it,norm_CTEF(it),Ekin_CTEF(it),Ebath_CTEF(it),Ecoup_CTEF(it)
        end do
        close(21)
      end if

    end subroutine CTEF_dynamics
!-----------------------------------------------------------------------------------------
    subroutine CTEF_dynamics_kernel
      implicit none
      complex(8),allocatable :: zK1_l(:,:,:,:,:),zK3_l(:,:,:,:,:)
      complex(8),allocatable :: zK1_sl(:,:,:,:,:),zK3_sl(:,:,:,:,:)
      complex(8),allocatable :: zK1_phase_ave(:,:,:,:,:),zK3_phase_ave(:,:,:,:,:)
      complex(8),allocatable :: zK1_t(:,:,:),zK3_t(:,:,:)
      complex(8) :: zpsi_store(Lsite,2),zHO_store(Lsite,2)
      complex(8) :: zweight, zweight0
      integer :: itraj_t, istore, itraj
      integer :: ntraj_tot_l, ntraj_tot
      integer :: ntraj_stable_l, ntraj_stable
      integer :: i_dm, j_dm, iphase
      real(8) :: norm
      real(8) :: phi0, phi
      integer,parameter :: ran_len = 1
      real(8) :: rvec(ran_len)
      logical :: is_norm_converged
      real(8) :: norm_CTEF_t(0:Nt+1)


      call allocate_GQME_kernel
      allocate(zK1_l(Lsite,Lsite,1,Lsite,0:Nt),zK3_l(Lsite,Lsite,1,Lsite,0:Nt))
      allocate(zK1_sl(Lsite,Lsite,1,Lsite,0:Nt),zK3_sl(Lsite,Lsite,1,Lsite,0:Nt))
      allocate(zK1_phase_ave(Lsite,Lsite,1,Lsite,0:Nt),zK3_phase_ave(Lsite,Lsite,1,Lsite,0:Nt))
      allocate(zK1_t(Lsite,Lsite,0:Nt),zK3_t(Lsite,Lsite,0:Nt))
      zK1_l = 0d0; zK3_l = 0d0

      ntraj_tot_l = 0; ntraj_stable_l = 0
      itraj = 0
      do itraj_t  = 1, max(Ntraj/nsize_store,1)
        zK1_sl = 0d0; zK3_sl = 0d0
        do istore = 1, min(nsize_store,Ntraj)
          itraj = itraj + 1
          call bath_sampling_correlated_gaussian(zHO_store,zweight0)

          call ranlux_double (rvec, ran_len)
          phi0 = rvec(1); phi0 = 2d0*pi*phi0
          if(myrank == 0 .and. mod(itraj,max(Ntraj/200,1))==0)write(*,*)"itraj=",itraj,"/",Ntraj
          if(mod(itraj,Nprocs) /= myrank)cycle

! == localized init wf
          i_dm = 1
          do j_dm = 1,Lsite
            zpsi_store = 0d0 
            zpsi_store(j_dm,1) = 1d0; zpsi_store(i_dm,2) = 1d0

            zweight = zweight0 * (conjg(zHO_store(i_dm,2))-zHO_store(j_dm,1))

            
            zK1_phase_ave = 0d0; zK3_phase_ave = 0d0
            is_norm_converged = .true.
            do iphase = 1,Nphase
              if(.not. is_norm_converged)exit

              phi = phi0 + 2d0*pi*dble(iphase-1)/Nphase
              call set_initial_condition(zpsi_store,zHO_store, &
                zpsi_CTEF, zHO_CTEF, phi, norm)

              call propagation_kernel(norm_CTEF_t,zK1_t, zK3_t)
              zK1_phase_ave(:,:,i_dm,j_dm,:) = zK1_phase_ave(:,:,i_dm,j_dm,:) &
                + gamma**2*zK1_t(:,:,:)*exp(-zI*phi)*norm*zweight
              zK3_phase_ave(:,:,i_dm,j_dm,:) = zK3_phase_ave(:,:,i_dm,j_dm,:) &
                - gamma * zK3_t(:,:,:)*exp(-zI*phi)*norm*zweight

              if(.not. abs(norm_CTEF_t(Nt)-1d0) < epsilon_norm )is_norm_converged = .false.

            end do

            if(is_norm_converged)then
              ntraj_stable_l = ntraj_stable_l + 1
              zK1_sl = zK1_sl + zK1_phase_ave
              zK3_sl = zK3_sl + zK3_phase_ave
            end if
            ntraj_tot_l = ntraj_tot_l + 1
            

          end do
        end do
        
        zK1_l = zK1_l + zK1_sl/dble(nsize_store)
        zK3_l = zK3_l + zK3_sl/dble(nsize_store)
      end do

      call MPI_ALLREDUCE(ntraj_tot_l,ntraj_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ntraj_stable_l,ntraj_stable,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myrank == 0)then
        write(*,*)"# of total   trajectries",ntraj_tot
        write(*,*)"# of stable  trajectries",ntraj_stable
        write(*,*)"# of skipped trajectries",ntraj_tot - ntraj_stable
      end if
      
      zK1_l = zK1_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)
      zK3_l = zK3_l&
        /dble(ntraj_stable)/dble(Nphase)*dble(Lsite)*dble(nsize_store)

      call MPI_ALLREDUCE(zK1_l,zK1,(Nt+1)*Lsite**3,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zK3_l,zK3,(Nt+1)*Lsite**3,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

      call refine_GQME_kernel_K1K3  ! Inversion symmetry, Unitarity, ... etc.
      call evaluate_GQME_kernel_full

      if(myrank == 0)then
        open(nfile_full_kernel,file=trim(file_full_kernel),form='unformatted')
        write(nfile_full_kernel)zK_full,zK1,zK3
        close(nfile_full_kernel)
      end if


    end subroutine CTEF_dynamics_kernel
!-----------------------------------------------------------------------------------------
    subroutine bath_sampling_correlated_gaussian(zHO_out,zweight)
      implicit none
      complex(8),intent(out) :: zHO_out(Lsite,2)
      complex(8),intent(out) :: zweight
      integer :: i
      real(8) :: x1,x2,p1,p2
      complex(8) :: z1,z2
      real(8) :: sigma, norm_fact, alpha, sigma1, sigma2

      sigma = sigma_correlated_gaussian
      norm_fact = 4d0*pi**2/( &
        (1d0+1d0/sigma**2) * &
        (1d0+1d0/(1d0+sigma**2)) )
      alpha = 1d0/(sigma**2 + 1d0)
      sigma1 = 1d0+alpha**2+(1d0-alpha)**2/sigma**2; sigma1 = 1d0/sqrt(sigma1)
      sigma2 = 1d0+1d0/sigma**2; sigma2 = 1d0/sqrt(sigma2)

      do i = 1,Lsite
        call gaussian_random_number(x1,p1)
        call gaussian_random_number(x2,p2)
        z1 = (x1 + zI * p1)*sigma1
        z2 = (x2 + zI * p2)*sigma2
        zHO_out(i,1) = z1
        zHO_out(i,2) = z2 + alpha*z1
      end do

      zweight = 1d0
      do i = 1,Lsite
        zweight = zweight * norm_fact/pi**2*exp( &
          0.5d0/sigma**2*abs(zHO_out(i,1)-zHO_out(i,2))**2 )
      end do

    end subroutine bath_sampling_correlated_gaussian
!-----------------------------------------------------------------------------------------
    subroutine bath_sampling_independent_gaussian(zHO_out,zweight)
      implicit none
      complex(8),intent(out) :: zHO_out(Lsite,2)
      complex(8),intent(out) :: zweight
      integer :: i
      real(8) :: x1,x2,p1,p2
      complex(8) :: z1,z2
      real(8) :: sigma

      sigma = sigma_independent_gaussian

      do i = 1,Lsite
        call gaussian_random_number(x1,p1)
        call gaussian_random_number(x2,p2)
        z1 = (x1 + zI * p1)*sqrt(0.5d0)
        z2 = (x2 + zI * p2)*sigma
        zHO_out(i,1) = z1
        zHO_out(i,2) = z2 + z1
      end do

      zweight = 1d0
      do i = 1,Lsite
        zweight = zweight * 2d0*sigma**2*exp( &
           0.5d0*abs(zHO_out(i,1))**2          &
          -0.5d0*abs(zHO_out(i,2))**2          &
          +0.5d0/sigma**2*abs(zHO_out(i,1)-zHO_out(i,2))**2 )
      end do

    end subroutine bath_sampling_independent_gaussian
!-----------------------------------------------------------------------------------------
    subroutine calc_zweight(zHO_in,zweight)
      implicit none
      complex(8),intent(in) :: zHO_in(Lsite,2)
      complex(8),intent(out) :: zweight
      integer :: i

      zweight = 1d0
      do i = 1, Lsite
        zweight = zweight * (4d0/3d0)*exp( &
          +0.5d0*abs(zHO_in(i,1) - zHO_in(i,2))**2 &
          )
      end do

    end subroutine calc_zweight
!-----------------------------------------------------------------------------------------
    subroutine set_initial_condition(zpsi_in,zHO_in,zpsi_out,zHO_out,phi,norm)
      implicit none
      complex(8),intent(in) :: zpsi_in(Lsite,2), zHO_in(Lsite,2)
      complex(8),intent(out) :: zpsi_out(Lsite,2), zHO_out(Lsite,2)
      real(8),intent(in) :: phi
      real(8),intent(out) :: norm
      complex(8) :: zs_elec, zs_bath

      zpsi_out(:,1) = zpsi_in(:,1)
      zpsi_out(:,2) = exp(zI*phi)*zpsi_in(:,2)
      zHO_out = zHO_in

      call calc_norm(zpsi_out,zHO_out,norm)

      zpsi_out = zpsi_out/sqrt(norm)

    end subroutine set_initial_condition
!-----------------------------------------------------------------------------------------
    subroutine propagation(norm_CTEF_t,Ekin_CTEF_t,Ebath_CTEF_t,Ecoup_CTEF_t)
      implicit none
      real(8),intent(out) :: norm_CTEF_t(0:Nt+1),Ekin_CTEF_t(0:Nt+1)
      real(8),intent(out) :: Ebath_CTEF_t(0:Nt+1),Ecoup_CTEF_t(0:Nt+1)
      integer :: it

      call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(0))

      zHO_dot_CTEF = 0d0
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      Ekin_CTEF_t(0) = sum(zEs_CTEF*zSb_CTEF)
      Ebath_CTEF_t(0) = sum(zEb_CTEF*zSs_CTEF)
      Ecoup_CTEF_t(0) = sum(zEc_CTEF)


      do it = 0,Nt-1

!        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call dt_evolve_Runge_Kutta(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(it+1))
        Ekin_CTEF_t(it+1) = sum(zEs_CTEF*zSb_CTEF)
        Ebath_CTEF_t(it+1) = sum(zEb_CTEF*zSs_CTEF)
        Ecoup_CTEF_t(it+1) = sum(zEc_CTEF)

        
      end do

    end subroutine propagation
!-----------------------------------------------------------------------------------------
    subroutine propagation_kernel(norm_CTEF_t,zK1_t,zK3_t)
      implicit none
      real(8),intent(out) :: norm_CTEF_t(0:Nt+1)
      complex(8),intent(out) :: zK1_t(Lsite,Lsite,0:Nt)
      complex(8),intent(out) :: zK3_t(Lsite,Lsite,0:Nt)
      complex(8) :: zrho_dm_t(Lsite,Lsite,2,2)
      integer :: it

      call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(0))

      zHO_dot_CTEF = 0d0
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
      call evaluate_kernel(0)

      do it = 0,Nt-1

!        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call dt_evolve_Runge_Kutta(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(it+1))
        call evaluate_kernel(it+1)
        
      end do
      
      contains
        subroutine evaluate_kernel(it_t)
          implicit none
          integer,intent(in) :: it_t
          integer :: i,j

          do i = 1, Lsite
            do j = 1,Lsite
              zrho_dm_t(i,j,1,1) = zpsi_CTEF(i,1)*conjg(zpsi_CTEF(j,1))*zSb_CTEF(1,1)
              zrho_dm_t(i,j,1,2) = zpsi_CTEF(i,2)*conjg(zpsi_CTEF(j,1))*zSb_CTEF(1,2)
              zrho_dm_t(i,j,2,1) = zpsi_CTEF(i,1)*conjg(zpsi_CTEF(j,2))*zSb_CTEF(2,1)
              zrho_dm_t(i,j,2,2) = zpsi_CTEF(i,2)*conjg(zpsi_CTEF(j,2))*zSb_CTEF(2,2)
            end do
          end do

          do i = 1,Lsite
            do j = 1,Lsite
              zK3_t(i,j,it_t) = sum(zrho_dm_t(i,j,:,:))
              zK1_t(i,j,it_t) = zrho_dm_t(i,j,1,1)*(&
                conjg(zHO_CTEF(i,1))+zHO_CTEF(i,1)-conjg(zHO_CTEF(j,1))-zHO_CTEF(j,1) ) &
                               +zrho_dm_t(i,j,2,2)*(&
                conjg(zHO_CTEF(i,2))+zHO_CTEF(i,2)-conjg(zHO_CTEF(j,2))-zHO_CTEF(j,2) ) &
                               +zrho_dm_t(i,j,1,2)*(&
                conjg(zHO_CTEF(i,1))+zHO_CTEF(i,2)-conjg(zHO_CTEF(j,1))-zHO_CTEF(j,2) ) &
                               +zrho_dm_t(i,j,2,1)*(&
                conjg(zHO_CTEF(i,2))+zHO_CTEF(i,1)-conjg(zHO_CTEF(j,2))-zHO_CTEF(j,1) ) 
            end do
          end do
        end subroutine evaluate_kernel

    end subroutine propagation_kernel
!-----------------------------------------------------------------------------------------
    subroutine calc_norm(zpsi_in,zHO_in,norm)
      implicit none
      complex(8),intent(in) :: zpsi_in(Lsite,2),zHO_in(Lsite,2)
      real(8),intent(out) :: norm
      complex(8) :: zs_elec, zs_bath

      zs_elec = sum( conjg(zpsi_in(:,1)) * zpsi_in(:,2) )
      zs_bath = sum( -0.5d0*abs(zHO_in(:,1))**2 -0.5d0*abs(zHO_in(:,2))**2 &
        +conjg(zHO_in(:,1))*zHO_in(:,2) )
      zs_bath = exp(zs_bath)

      norm = sum(abs(zpsi_in(:,:))**2) + 2d0*real(zs_elec*zs_bath)

    end subroutine calc_norm
!-----------------------------------------------------------------------------------------
    subroutine refine_effective_hamiltonian(zpsi_in,zHO_in,zHO_dot_inout)
      implicit none
      complex(8),intent(in) :: zpsi_in(Lsite,2),zHO_in(Lsite,2)
      complex(8),intent(inout) :: zHO_dot_inout(Lsite,2)
      complex(8) :: zhpsi_t(Lsite,2)
      complex(8) :: zs, zvec(2)
      integer :: i,j, iscf

! zSs
      zSs_CTEF(1,1) = sum(abs(zpsi_in(:,1))**2)
      zSs_CTEF(2,2) = sum(abs(zpsi_in(:,2))**2)
      zSs_CTEF(1,2) = sum(conjg(zpsi_in(:,1))*zpsi_in(:,2))
      zSs_CTEF(2,1) = conjg(zSs_CTEF(1,2))

! zSb
      zSb_CTEF(1,1) = 1d0; zSb_CTEF(2,2) = 1d0
      zs = -0.5d0*sum(abs(zHO_in(:,:))**2)+sum(conjg(zHO_in(:,1))*zHO_in(:,2))
      zSb_CTEF(1,2) = exp(zs)
      zSb_CTEF(2,1) = conjg(zSb_CTEF(1,2))
      call inverse_2x2_matrix(zSb_CTEF,zSb_inv_CTEF)

      zSsb_CTEF = zSs_CTEF*zSb_CTEF

      call inverse_2x2_matrix(zSsb_CTEF,zSsb_inv_CTEF)

! zX_HO
      do i = 1,2
        do j = 1,2
          zX_HO_CTEF(:,i,j) = (conjg(zHO_in(:,i)) + zHO_in(:,j)) &
            *zSb_CTEF(i,j)
        end do
      end do

! zEb
      do i = 1,2
        do j = 1,2
          zEb_CTEF(i,j) = omega0*( &
            sum( conjg(zHO_in(:,i))*zHO_in(:,j) ) + 0.5d0*dble(Lsite) &
            )*zSb_CTEF(i,j)
        end do
      end do

! zEs
      call hs_zpsi(zpsi_in,zhpsi_t)
      zEs_CTEF(1,1) = real( sum( conjg(zpsi_in(:,1))*zhpsi_t(:,1) ))
      zEs_CTEF(2,2) = real( sum( conjg(zpsi_in(:,2))*zhpsi_t(:,2) ))
      zEs_CTEF(1,2) = sum(conjg(zpsi_in(:,1))*zhpsi_t(:,2))
      zEs_CTEF(2,1) = conjg(zEs_CTEF(1,2))

! zEc
      zEc_CTEF(1,1) = sum( zX_HO_CTEF(:,1,1)*abs(zpsi_in(:,1))**2)
      zEc_CTEF(2,2) = sum( zX_HO_CTEF(:,2,2)*abs(zpsi_in(:,2))**2)
      zEc_CTEF(1,2) = sum( zX_HO_CTEF(:,1,2)*conjg(zpsi_in(:,1))*zpsi_in(:,2))
      zEc_CTEF(2,1) = conjg(zEc_CTEF(1,2))
      zEc_CTEF = -gamma* zEc_CTEF 

      do i = 1, Lsite
        zF_HO_CTEF(i,1) = abs(zpsi_in(i,1))**2 &
          + conjg(zpsi_in(i,1))*zpsi_in(i,2)*zSb_CTEF(1,2)
        zF_HO_CTEF(i,2) = abs(zpsi_in(i,2))**2 &
          + conjg(zpsi_in(i,2))*zpsi_in(i,1)*zSb_CTEF(2,1)
      end do
      zF_HO_CTEF = -gamma* zF_HO_CTEF


      do iscf = 1, Nscf_refine

! zDb
        zDb_CTEF(1,1) = real( 0.5d0*zI*sum(conjg(zHO_in(:,1))*zHO_dot_inout(:,1) &
          -conjg(zHO_dot_inout(:,1))*zHO_in(:,1) ) )
        zDb_CTEF(2,2) = real( 0.5d0*zI*sum(conjg(zHO_in(:,2))*zHO_dot_inout(:,2) &
          -conjg(zHO_dot_inout(:,2))*zHO_in(:,2) ) )
        zDb_CTEF(1,2) = zI*sum(-0.5d0*( &
          conjg(zHO_dot_inout(:,2))*zHO_in(:,2) &
          +zHO_dot_inout(:,2)*conjg(zHO_in(:,2)) ) &
          +conjg(zHO_in(:,1))*zHO_dot_inout(:,2) )*zSb_CTEF(1,2)
        zDb_CTEF(2,1) = zI*sum(-0.5d0*( &
          conjg(zHO_dot_inout(:,1))*zHO_in(:,1) &
          +zHO_dot_inout(:,1)*conjg(zHO_in(:,1)) ) &
          +conjg(zHO_in(:,2))*zHO_dot_inout(:,1) )*zSb_CTEF(2,1)

        call heff_zpsi(zpsi_in,zhpsi_t)
        zDs_CTEF(1,1) = sum( conjg(zpsi_in(:,1))*zhpsi_t(:,1) )
        zDs_CTEF(2,2) = sum( conjg(zpsi_in(:,2))*zhpsi_t(:,2) )
        zDs_CTEF(1,2) = sum( conjg(zpsi_in(:,1))*zhpsi_t(:,2) )
        zDs_CTEF(2,1) = sum( conjg(zpsi_in(:,2))*zhpsi_t(:,1) )


        zHb_eff_CTEF = omega0*zSsb_CTEF
        zHb_eff_CTEF(1,1) = zHb_eff_CTEF(1,1) - zI*real(zDs_CTEF(1,1)/zI) &
          +real(zDs_CTEF(1,2)*zSb_CTEF(1,2) + zSs_CTEF(1,2)*zDb_CTEF(1,2)) &
          -real(zEs_CTEF(1,2)*zSb_CTEF(1,2) + zEb_CTEF(1,2)*zSs_CTEF(1,2) + zEc_CTEF(1,2))

        zHb_eff_CTEF(2,2) = zHb_eff_CTEF(2,2) - zI*real(zDs_CTEF(2,2)/zI) &
          +real(zDs_CTEF(2,1)*zSb_CTEF(2,1) + zSs_CTEF(2,1)*zDb_CTEF(2,1)) &
          -real(zEs_CTEF(2,1)*zSb_CTEF(2,1) + zEb_CTEF(2,1)*zSs_CTEF(2,1) + zEc_CTEF(2,1))

        zHb_eff_CTEF(1,2) = zHb_eff_CTEF(1,2) &
          -zDs_CTEF(1,2)*zSb_CTEF(1,2) -zSs_CTEF(1,2)*zDb_CTEF(1,2) &
          +zEs_CTEF(1,2)*zSb_CTEF(1,2) + zEb_CTEF(1,2)*zSs_CTEF(1,2) + zEc_CTEF(1,2)

        zHb_eff_CTEF(2,1) = zHb_eff_CTEF(2,1) &
          -zDs_CTEF(2,1)*zSb_CTEF(2,1) -zSs_CTEF(2,1)*zDb_CTEF(2,1) &
          +zEs_CTEF(2,1)*zSb_CTEF(2,1) + zEb_CTEF(2,1)*zSs_CTEF(2,1) + zEc_CTEF(2,1)


        do i = 1, Lsite
          zvec(1) = zHb_eff_CTEF(1,1)*zHO_in(i,1) + zHb_eff_CTEF(1,2)*zHO_in(i,2) &
            + zF_HO_CTEF(i,1)
          zvec(2) = zHb_eff_CTEF(2,1)*zHO_in(i,1) + zHb_eff_CTEF(2,2)*zHO_in(i,2) &
            + zF_HO_CTEF(i,2)
          zHO_dot_inout(i,1) = (zSsb_inv_CTEF(1,1)*zvec(1) + zSsb_inv_CTEF(1,2)*zvec(2) )/zI 
          zHO_dot_inout(i,2) = (zSsb_inv_CTEF(2,1)*zvec(1) + zSsb_inv_CTEF(2,2)*zvec(2) )/zI 
        end do

      end do

    end subroutine refine_effective_hamiltonian
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_Runge_Kutta(zpsi_inout,zHO_inout,zHO_dot_inout)
      implicit none
      complex(8),intent(inout) :: zpsi_inout(Lsite,2),zHO_inout(Lsite,2)
      complex(8),intent(inout) :: zHO_dot_inout(Lsite,2)
      complex(8) :: zpsi_t0(Lsite,2),zHO_t0(Lsite,2)
      complex(8) :: zpsi_t1(Lsite,2),zHO_t1(Lsite,2)
      complex(8) :: zpsi_t2(Lsite,2),zHO_t2(Lsite,2)
      complex(8) :: zpsi_t3(Lsite,2),zHO_t3(Lsite,2)
      complex(8) :: zpsi_t4(Lsite,2),zHO_t4(Lsite,2)
      integer :: iscf

      zpsi_t0 = zpsi_inout
      zHO_t0 = zHO_inout

! k1
      call heff_zpsi(zpsi_inout,zpsi_t1)
      zpsi_t1 = -zI*zpsi_t1
      zHO_t1  = zHO_dot_inout

! k2
      zpsi_inout = zpsi_inout + 0.5d0*dt*zpsi_t1
      zHO_inout  = zHO_inout  + 0.5d0*dt*zHO_t1
      call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
      call heff_zpsi(zpsi_inout,zpsi_t2)
      zpsi_t2 = -zI*zpsi_t2
      zHO_t2  = zHO_dot_inout

! k3
      zpsi_inout = zpsi_t0 + 0.5d0*dt*zpsi_t2
      zHO_inout  = zHO_t0  + 0.5d0*dt*zHO_t2
      call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
      call heff_zpsi(zpsi_inout,zpsi_t3)
      zpsi_t3 = -zI*zpsi_t3
      zHO_t3  = zHO_dot_inout

! k4
      zpsi_inout = zpsi_t0 + dt*zpsi_t3
      zHO_inout  = zHO_t0  + dt*zHO_t3
      call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
      call heff_zpsi(zpsi_inout,zpsi_t4)
      zpsi_t4 = -zI*zpsi_t4
      zHO_t4  = zHO_dot_inout

      zpsi_inout = zpsi_t0 + dt/6d0*(zpsi_t1 + 2d0*zpsi_t2 + 2d0*zpsi_t3 + zpsi_t4)
      zHO_inout  = zHO_t0  + dt/6d0*(zHO_t1  + 2d0*zHO_t2  + 2d0*zHO_t3  + zHO_t4)

      call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)

    end subroutine dt_evolve_Runge_Kutta
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_etrs(zpsi_inout,zHO_inout,zHO_dot_inout)
      implicit none
      complex(8),intent(inout) :: zpsi_inout(Lsite,2),zHO_inout(Lsite,2)
      complex(8),intent(inout) :: zHO_dot_inout(Lsite,2)
      complex(8) :: zpsi_t(Lsite,2),zHO_t(Lsite,2)
      integer :: iscf

! t -> t + dt/2
      call dt_evolve_elec(zpsi_inout,dt*0.5d0)
      call dt_evolve_bath(zHO_inout,zHO_dot_inout,dt*0.5d0)
!      call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)
      zpsi_t = zpsi_inout
      zHO_t = zHO_inout


      call dt_evolve_elec(zpsi_inout,dt*0.5d0)
      call dt_evolve_bath(zHO_inout,zHO_dot_inout,dt*0.5d0)
!      call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)

      do iscf = 1, Nscf_pred_corr
        call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
        zHO_inout = zHO_t
        zpsi_inout = zpsi_t

        call dt_evolve_elec(zpsi_inout,dt*0.5d0)
        call dt_evolve_bath(zHO_inout,zHO_dot_inout,dt*0.5d0)
!        call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)

      end do

    end subroutine dt_evolve_etrs
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_elec(zpsi_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zpsi_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Nexp_Taylor = 6
      complex(8) :: zpsi_t(Lsite,2),zhpsi_t(Lsite,2)
      integer :: iexp
      complex(8) :: zfact

      zpsi_t = zpsi_inout
      zfact = 1d0
      do iexp = 1,Nexp_Taylor
        zfact = zfact*(-zI*dt_t)/iexp
        call heff_zpsi(zpsi_t,zhpsi_t)
        zpsi_inout = zpsi_inout + zfact*zhpsi_t
        zpsi_t = zhpsi_t
      end do

    end subroutine dt_evolve_elec
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath(zHO_inout,zHO_dot_in,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      complex(8),intent(in) :: zHO_dot_in(Lsite,2)

      select case(iflag_bath_propagator)
      case(bath_propagator_direct)
        call dt_evolve_bath_direct(zHO_inout,zHO_dot_in,dt_t)
      case(bath_propagator_taylor)
        call dt_evolve_bath_taylor(zHO_inout,dt_t)
      case(bath_propagator_taylor_mod)
        call dt_evolve_bath_taylor_mod(zHO_inout,dt_t)
      case(bath_propagator_diag)
        call dt_evolve_bath_diag(zHO_inout,dt_t)
      case(bath_propagator_CrankNicolson)
        call dt_evolve_bath_CrankNicolson(zHO_inout,dt_t)
      case default
        stop 'wrong bath propagator'
      end select

    end subroutine dt_evolve_bath
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath_direct(zHO_inout,zHO_dot_in,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      complex(8),intent(in) :: zHO_dot_in(Lsite,2)
      real(8),intent(in) :: dt_t

      zHO_inout = zHO_inout + dt_t*zHO_dot_CTEF
      return
    end subroutine dt_evolve_bath_direct
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath_taylor(zHO_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Nexp_Taylor = 6
      complex(8) :: zHO_t(Lsite,2),zhHO_t(Lsite,2), zF_HO_eff(Lsite,2)
      complex(8) :: zHeff(2,2)
      integer :: iexp
      complex(8) :: zfact

      zF_HO_eff(:,1) = zSsb_inv_CTEF(1,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(1,2)*zF_HO_CTEF(:,2)
      zF_HO_eff(:,2) = zSsb_inv_CTEF(2,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(2,2)*zF_HO_CTEF(:,2)
      zHeff = matmul(zSsb_inv_CTEF,zHb_eff_CTEF)

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

      zHO_t = zHO_inout
      zfact = 1d0
      do iexp = 1,Nexp_Taylor
        zfact = zfact*(-zI*dt_t)/iexp
        zhHO_t(:,1) = zHeff(1,1)*zHO_t(:,1) + zHeff(1,2)*zHO_t(:,2)
        zhHO_t(:,2) = zHeff(2,1)*zHO_t(:,1) + zHeff(2,2)*zHO_t(:,2)

        zHO_inout = zHO_inout + zfact*zhHO_t
        zHO_t = zhHO_t
      end do

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

    end subroutine dt_evolve_bath_taylor
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath_taylor_mod(zHO_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Nexp_Taylor = 6 !, Nexp_mult = 0
      complex(8) :: zHO_t(Lsite,2),zhHO_t(Lsite,2), zF_HO_eff(Lsite,2)
      complex(8) :: zHeff(2,2)
      complex(8) :: zvec_t(Lsite,2)
      complex(8) :: zexp_Hmat(2,2),zmat1(2,2),zmat2(2,2)
      complex(8) :: zs
      integer :: iexp
      real(8) :: ss

      zF_HO_eff(:,1) = zSsb_inv_CTEF(1,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(1,2)*zF_HO_CTEF(:,2)
      zF_HO_eff(:,2) = zSsb_inv_CTEF(2,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(2,2)*zF_HO_CTEF(:,2)
      zHeff = matmul(zSsb_inv_CTEF,zHb_eff_CTEF)
      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

!      zmat1 = -zI*dt_t*0.5d0**Nexp_mult*zHeff
      zmat1 = -zI*dt_t*zHeff
      zexp_Hmat(:,1) = (/1d0,0d0/)
      zexp_Hmat(:,2) = (/0d0,1d0/)
      zexp_Hmat = zexp_Hmat + zmat1
      zmat2 = zmat1
      ss = 1d0

      do iexp = 2,Nexp_Taylor
        ss = ss/iexp
        zmat2 = matmul(zmat2,zmat1)
        zexp_Hmat = zexp_Hmat + ss * zmat2
      end do

!      do iexp = 1,Nexp_mult
!        zexp_Hmat = matmul(zexp_Hmat,zexp_Hmat)
!      end do

      
      zvec_t = zHO_inout

      zHO_inout(:,1) = zexp_Hmat(1,1)*zvec_t(:,1) + zexp_Hmat(1,2)*zvec_t(:,2)
      zHO_inout(:,2) = zexp_Hmat(2,1)*zvec_t(:,1) + zexp_Hmat(2,2)*zvec_t(:,2)

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

    end subroutine dt_evolve_bath_taylor_mod
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath_diag(zHO_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      real(8),parameter :: epsilon = 1d-6
      integer,parameter :: Nexp_Taylor = 6
      complex(8) :: zlambda(2), zeig_vec(2,2), zeig_vec_inv(2,2)
      complex(8) :: zHO_t(Lsite,2),zhHO_t(Lsite,2), zF_HO_eff(Lsite,2)
      complex(8) :: zHeff(2,2)
      complex(8) :: zvec_t(Lsite,2),zexp_Hmat(2,2),zmat(2,2),zmat2(2,2)
      integer :: iexp
      complex(8) :: zs
      real(8) :: ss

      zF_HO_eff(:,1) = zSsb_inv_CTEF(1,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(1,2)*zF_HO_CTEF(:,2)
      zF_HO_eff(:,2) = zSsb_inv_CTEF(2,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(2,2)*zF_HO_CTEF(:,2)
      zHeff = matmul(zSsb_inv_CTEF,zHb_eff_CTEF)

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

      zs = (zHeff(1,1) - zHeff(2,2))**2 &
        + 4d0*zHeff(1,2)*zHeff(2,1)
      zlambda(1) = 0.5d0*( zHeff(1,1) + zHeff(2,2) + sqrt(zs) )
      zlambda(2) = 0.5d0*( zHeff(1,1) + zHeff(2,2) - sqrt(zs) )

      if(abs(zHeff(1,2)*zHeff(2,1)) > epsilon &
        .and. abs(zlambda(1)*zlambda(2)) >epsilon)then

        if( abs(zHeff(1,1) -zlambda(1)) < abs(zHeff(2,2) -zlambda(1)) )then
          zeig_vec(1,1) = 1d0
          zeig_vec(2,1) = (zlambda(1) - zHeff(1,1))/zHeff(1,2)
          zeig_vec(2,2) = 1d0
          zeig_vec(1,2) = (zlambda(2) - zHeff(2,2))/zHeff(2,1)
        else
          zeig_vec(1,2) = 1d0
          zeig_vec(2,2) = (zlambda(2) - zHeff(1,1))/zHeff(1,2)
          zeig_vec(2,1) = 1d0
          zeig_vec(1,1) = (zlambda(1) - zHeff(2,2))/zHeff(2,1)
        end if

        call inverse_2x2_matrix(zeig_vec,zeig_vec_inv)
        zexp_Hmat(1,1) = exp(-zI*dt_t*zlambda(1))
        zexp_Hmat(2,1) = 0d0
        zexp_Hmat(1,2) = 0d0
        zexp_Hmat(2,2) = exp(-zI*dt_t*zlambda(2))
        zmat = matmul(zexp_Hmat,zeig_vec_inv)
        zexp_Hmat = matmul(zeig_vec,zmat)

      else
        zmat(1,1) = 0d0; zmat(2,1) = -zI*dt_t*zHeff(2,1)
        zmat(1,2) = -zI*dt_t*zHeff(1,2); zmat(2,2) = 0d0
        zexp_Hmat(1,1) = 1d0
        zexp_Hmat(2,1) = 0d0
        zexp_Hmat(1,2) = 0d0
        zexp_Hmat(2,2) = 1d0
        zexp_Hmat = zexp_Hmat +zmat
        zmat2 = zmat
        ss = 1d0
        do iexp = 2,Nexp_Taylor
          ss = ss/iexp
          zmat2 = matmul(zmat2,zmat)
          zexp_Hmat = zexp_Hmat + ss*zmat2
        end do

        zmat(1,1) = exp(-zI*0.5d0*dt_t*zHeff(1,1)); zmat(2,1) = 0d0
        zmat(1,2) = 0d0; zmat(2,2) = exp(-zI*0.5d0*dt_t*zHeff(2,2))
        zexp_Hmat = matmul(zexp_Hmat,zmat)
        zexp_Hmat = matmul(zmat,zexp_Hmat)

      end if

      
      zvec_t = zHO_inout

      zHO_inout(:,1) = zexp_Hmat(1,1)*zvec_t(:,1) + zexp_Hmat(1,2)*zvec_t(:,2)
      zHO_inout(:,2) = zexp_Hmat(2,1)*zvec_t(:,1) + zexp_Hmat(2,2)*zvec_t(:,2)

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

    end subroutine dt_evolve_bath_diag
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath_CrankNicolson(zHO_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      real(8),parameter :: umat(2,2) = reshape( (/1d0, 0d0, 0d0, 1d0/), (/2,2/) )
      complex(8) :: zmat1(2,2),zmat2(2,2),zmat2_inv(2,2)
      complex(8) :: zvec_t(Lsite,2)
      complex(8) :: zHeff(2,2),zF_HO_eff(Lsite,2)

      zF_HO_eff(:,1) = zSsb_inv_CTEF(1,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(1,2)*zF_HO_CTEF(:,2)
      zF_HO_eff(:,2) = zSsb_inv_CTEF(2,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(2,2)*zF_HO_CTEF(:,2)
      zHeff  = matmul(zSsb_inv_CTEF,zHb_eff_CTEF)
      zmat1 = umat - 0.5d0*zi*dt_t*zHeff
      zmat2 = umat + 0.5d0*zi*dt_t*zHeff
      call inverse_2x2_matrix(zmat2,zmat2_inv)

      zvec_t(:,1) = zmat1(1,1)*zHO_inout(:,1) + zmat1(1,2)*zHO_inout(:,2) &
        -zI*dt_t* zF_HO_eff(:,1)
      zvec_t(:,2) = zmat1(2,1)*zHO_inout(:,1) + zmat1(2,2)*zHO_inout(:,2) &
        -zI*dt_t* zF_HO_eff(:,2)

      zmat1 = (zI*zSsb_CTEF + 0.5d0*dt_t*zHb_eff_CTEF)
      zmat2 = (zI*zSsb_CTEF - 0.5d0*dt_t*zHb_eff_CTEF)


      zHO_inout(:,1) = zmat2_inv(1,1)*zvec_t(:,1) + zmat2_inv(1,2)*zvec_t(:,2)
      zHO_inout(:,2) = zmat2_inv(2,1)*zvec_t(:,1) + zmat2_inv(2,2)*zvec_t(:,2)

    end subroutine dt_evolve_bath_CrankNicolson
!-----------------------------------------------------------------------------------------
    subroutine hs_zpsi(zpsi_in,zhpsi_out)
      implicit none
      complex(8),intent(in) :: zpsi_in(Lsite,2)
      complex(8),intent(out) :: zhpsi_out(Lsite,2)
      integer :: i

      zhpsi_out(1,1) = -t0*(zpsi_in(2,1) + zpsi_in(Lsite,1))
      do i = 2,Lsite-1
        zhpsi_out(i,1) = -t0*(zpsi_in(i+1,1) + zpsi_in(i-1,1))
      end do
      zhpsi_out(Lsite,1) = -t0*(zpsi_in(1,1) + zpsi_in(Lsite-1,1))

      zhpsi_out(1,2) = -t0*(zpsi_in(2,2) + zpsi_in(Lsite,2))
      do i = 2,Lsite-1
        zhpsi_out(i,2) = -t0*(zpsi_in(i+1,2) + zpsi_in(i-1,2))
      end do
      zhpsi_out(Lsite,2) = -t0*(zpsi_in(1,2) + zpsi_in(Lsite-1,2))


    end subroutine hs_zpsi
!-----------------------------------------------------------------------------------------
    subroutine heff_zpsi(zpsi_in,zhpsi_out)
      implicit none
      complex(8),intent(in) :: zpsi_in(Lsite,2)
      complex(8),intent(out) :: zhpsi_out(Lsite,2)
      complex(8) :: zhpsi_t(Lsite,2)
      complex(8) :: zhs_psi_t(Lsite,2)

      call hs_zpsi(zpsi_in,zhs_psi_t)

      zhpsi_t(:,1) = zSb_CTEF(1,1)*zhs_psi_t(:,1) + zSb_CTEF(1,2)*zhs_psi_t(:,2)
      zhpsi_t(:,2) = zSb_CTEF(2,1)*zhs_psi_t(:,1) + zSb_CTEF(2,2)*zhs_psi_t(:,2)

      zhpsi_t(:,1) = zhpsi_t(:,1) &
        -gamma*zX_HO_CTEF(:,1,1)*zpsi_in(:,1) -gamma*zX_HO_CTEF(:,1,2)*zpsi_in(:,2)
      zhpsi_t(:,2) = zhpsi_t(:,2) &
        -gamma*zX_HO_CTEF(:,2,1)*zpsi_in(:,1) -gamma*zX_HO_CTEF(:,2,2)*zpsi_in(:,2)

      zhpsi_t(:,1) = zhpsi_t(:,1) &
        +zEb_CTEF(1,1)*zpsi_in(:,1) +zEb_CTEF(1,2)*zpsi_in(:,2)
      zhpsi_t(:,2) = zhpsi_t(:,2) &
        +zEb_CTEF(2,1)*zpsi_in(:,1) +zEb_CTEF(2,2)*zpsi_in(:,2)

      zhpsi_t(:,1) = zhpsi_t(:,1) &
        -zDb_CTEF(1,1)*zpsi_in(:,1) -zDb_CTEF(1,2)*zpsi_in(:,2)
      zhpsi_t(:,2) = zhpsi_t(:,2) &
        -zDb_CTEF(2,1)*zpsi_in(:,1) -zDb_CTEF(2,2)*zpsi_in(:,2)

      zhpsi_out(:,1) = zSb_inv_CTEF(1,1)*zhpsi_t(:,1) + zSb_inv_CTEF(1,2)*zhpsi_t(:,2)
      zhpsi_out(:,2) = zSb_inv_CTEF(2,1)*zhpsi_t(:,1) + zSb_inv_CTEF(2,2)*zhpsi_t(:,2)

    end subroutine heff_zpsi
!-----------------------------------------------------------------------------------------
end module CTEF_mod
