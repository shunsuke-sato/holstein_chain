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
  complex(8),allocatable :: zpsi_CTEF(:,:),zHO_CTEF(:,:)
  complex(8),allocatable :: zHO_dot_CTEF(:,:)

  complex(8) :: zSs_CTEF(2,2), zDs_CTEF(2,2)
  complex(8) :: zSb_inv_CTEF(2,2)
  complex(8) :: zSb_CTEF(2,2), zDb_CTEF(2,2)
  complex(8) :: zEs_CTEF(2,2), zEc_CTEF(2,2), zEb_CTEF(2,2)
  complex(8) :: zHb_eff_CTEF(2,2), zSsb_CTEF(2,2), zSsb_inv_CTEF(2,2)
  complex(8),allocatable :: zX_HO_CTEF(:,:,:), zF_HO_CTEF(:,:)

  public :: CTEF

  contains

!-----------------------------------------------------------------------------------------
    subroutine CTEF
      implicit none

      call CTEF_allocation
      call CTEF_dynamics

    end subroutine CTEF
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
      integer :: itraj, iphase, it, i,j
      integer :: i_antithetic, j_antithetic
      integer :: i_dm, j_dm
      real(8) :: norm_CTEF(0:Nt+1),norm_CTEF_l(0:Nt+1),norm_CTEF_t(0:Nt+1)
      real(8) :: norm_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ekin_CTEF(0:Nt+1),Ekin_CTEF_l(0:Nt+1),Ekin_CTEF_t(0:Nt+1)
      real(8) :: Ekin_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ebath_CTEF(0:Nt+1),Ebath_CTEF_l(0:Nt+1),Ebath_CTEF_t(0:Nt+1)
      real(8) :: Ebath_CTEF_phase_ave(0:Nt+1)
      real(8) :: Ecoup_CTEF(0:Nt+1),Ecoup_CTEF_l(0:Nt+1),Ecoup_CTEF_t(0:Nt+1)
      real(8) :: Ecoup_CTEF_phase_ave(0:Nt+1)
      complex(8) :: zrho_dm
      real(8) :: x1,x2,p1,p2
      integer,parameter :: ran_len = 1
      real(8) :: rvec(ran_len)

      norm_CTEF_l = 0d0; Ekin_CTEF_l = 0d0
      Ebath_CTEF_l = 0d0; Ecoup_CTEF_l = 0d0


      do itraj = 1, Ntraj

! == bath distribution
        do i = 1,Lsite
          call gaussian_random_number(x1,p1)
          call gaussian_random_number(x2,p2)
          zHO_gauss_store(i,1) = (x1 + zI * p1)*sqrt(2d0/3d0)
          zHO_gauss_store(i,2) = (x2 + zI * p2)*sqrt(0.5d0)
        end do
! == bath distribution


        CALL ranlux_double (rvec, ran_len)
        phi0 = rvec(1); phi0 = 2d0*pi*phi0
        if(myrank == 0 .and. mod(itraj,Ntraj/200)==0)write(*,*)"itraj=",itraj,"/",Ntraj
        if(mod(itraj,Nprocs) /= myrank)cycle
!        write(*,*)"itraj=",itraj,"/",Ntraj
        
        do i_antithetic = 1,4
          do j_antithetic = 1,4

            select case(i_antithetic)
            case(1)
              zHO_store(:,1) = zHO_gauss_store(:,1)
            case(2)
              zHO_store(:,1) = -zHO_gauss_store(:,1)
            case(3)
              zHO_store(:,1) = conjg(zHO_gauss_store(:,1))
            case(4)
              zHO_store(:,1) = -conjg(zHO_gauss_store(:,1))
            end select

            select case(j_antithetic)
            case(1)
              zHO_store(:,2) = zHO_gauss_store(:,2)+0.5d0*zHO_store(:,1)
            case(2)
              zHO_store(:,2) = -zHO_gauss_store(:,2)+0.5d0*zHO_store(:,1)
            case(3)
              zHO_store(:,2) = conjg(zHO_gauss_store(:,2))+0.5d0*zHO_store(:,1)
            case(4)
              zHO_store(:,2) = -conjg(zHO_gauss_store(:,2))+0.5d0*zHO_store(:,1)
            end select

            call calc_zweight(zHO_store,zweight0)

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
              do iphase = 1,Nphase
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


              end do
              norm_CTEF_l = norm_CTEF_l + norm_CTEF_phase_ave
              Ekin_CTEF_l = Ekin_CTEF_l + Ekin_CTEF_phase_ave
              Ebath_CTEF_l = Ebath_CTEF_l + Ebath_CTEF_phase_ave
              Ecoup_CTEF_l = Ecoup_CTEF_l + Ecoup_CTEF_phase_ave
              
            end do
          end do
        end do
        
      end do

      norm_CTEF_l = norm_CTEF_l/dble(Ntraj*Nphase)/16d0
      Ekin_CTEF_l = Ekin_CTEF_l/dble(Ntraj*Nphase)/16d0
      Ebath_CTEF_l = Ebath_CTEF_l/dble(Ntraj*Nphase)/16d0
      Ecoup_CTEF_l = Ecoup_CTEF_l/dble(Ntraj*Nphase)/16d0
      call MPI_ALLREDUCE(norm_CTEF_l,norm_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ekin_CTEF_l,Ekin_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ebath_CTEF_l,Ebath_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ecoup_CTEF_l,Ecoup_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myrank == 0)then
        open(21,file="CTEF_norm.out")
        do it = 0,Nt+1
          write(21,"(999e26.16e3)")dt*it,norm_CTEF(it),Ekin_CTEF(it),Ebath_CTEF(it),Ecoup_CTEF(it)
        end do
        close(21)
      end if

    end subroutine CTEF_dynamics
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
      Ekin_CTEF_t(0) = sum(zEs_CTEF*zSb_CTEF)
      Ebath_CTEF_t(0) = sum(zEb_CTEF*zSs_CTEF)
      Ecoup_CTEF_t(0) = sum(zEc_CTEF)


      do it = 0,Nt-1

        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(it+1))
        Ekin_CTEF_t(it+1) = sum(zEs_CTEF*zSb_CTEF)
        Ebath_CTEF_t(it+1) = sum(zEb_CTEF*zSs_CTEF)
        Ecoup_CTEF_t(it+1) = sum(zEc_CTEF)

        
      end do

    end subroutine propagation
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
    subroutine dt_evolve_etrs(zpsi_inout,zHO_inout,zHO_dot_inout)
      implicit none
      complex(8),intent(inout) :: zpsi_inout(Lsite,2),zHO_inout(Lsite,2)
      complex(8),intent(inout) :: zHO_dot_inout(Lsite,2)
      complex(8) :: zpsi_t(Lsite,2),zHO_t(Lsite,2)
      integer :: iscf

! t -> t + dt/2
      call dt_evolve_elec(zpsi_inout,dt*0.5d0)
!      call dt_evolve_bath_direct(zHO_inout,zHO_dot_inout,dt*0.5d0)
      call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)
      zpsi_t = zpsi_inout
      zHO_t = zHO_inout


      call dt_evolve_elec(zpsi_inout,dt*0.5d0)
!      call dt_evolve_bath_direct(zHO_inout,zHO_dot_inout,dt*0.5d0)
      call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)

      do iscf = 1, Nscf_pred_corr
        call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
        zHO_inout = zHO_t
        zpsi_inout = zpsi_t

        call dt_evolve_elec(zpsi_inout,dt*0.5d0)
!        call dt_evolve_bath_direct(zHO_inout,zHO_dot_inout,dt*0.5d0)
        call dt_evolve_bath_taylor(zHO_inout,dt*0.5d0)

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
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Nexp_Taylor = 6
      complex(8) :: zHO_t(Lsite,2),zhHO_t(Lsite,2), zF_HO_eff(Lsite,2)
      integer :: iexp
      complex(8) :: zfact

      zF_HO_eff(:,1) = zSsb_inv_CTEF(1,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(1,2)*zF_HO_CTEF(:,2)
      zF_HO_eff(:,2) = zSsb_inv_CTEF(2,1)*zF_HO_CTEF(:,1) &
                     + zSsb_inv_CTEF(2,2)*zF_HO_CTEF(:,2)

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

      zHO_t = zHO_inout
      zfact = 1d0
      do iexp = 1,Nexp_Taylor
        zfact = zfact*(-zI*dt_t)/iexp
        zhHO_t(:,1) = zHb_eff_CTEF(1,1)*zHO_t(:,1) + zHb_eff_CTEF(1,2)*zHO_t(:,2)
        zhHO_t(:,2) = zHb_eff_CTEF(2,1)*zHO_t(:,1) + zHb_eff_CTEF(2,2)*zHO_t(:,2)

        zHO_inout = zHO_inout + zfact*zhHO_t
        zHO_t = zhHO_t
      end do

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_eff

    end subroutine dt_evolve_bath_taylor
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
