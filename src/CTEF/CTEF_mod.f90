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
      complex(8) :: zpsi_store(Lsite,2),zHO_store(Lsite,2)
      complex(8) :: zweight
      real(8) :: norm, phi0,phi
      integer :: itraj, iphase, it
      real(8) :: norm_CTEF(0:Nt+1), Ekin_CTEF(0:Nt+1)
      real(8) :: norm_CTEF_l(0:Nt+1), Ekin_CTEF_l(0:Nt+1)
      real(8) :: norm_CTEF_t(0:Nt+1), Ekin_CTEF_t(0:Nt+1)

      norm_CTEF_l = 0d0; Ekin_CTEF = 0d0


      do itraj = 1, Ntraj

        call init_forward_backward_trajectries(zpsi_store,zHO_store,zweight)
        call random_number(phi0); phi0 = 2d0*pi*phi0
        if(mod(itraj,Nprocs) /= myrank)cycle
        if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj
!        write(*,*)"itraj=",itraj,"/",Ntraj

        do iphase = 1,Nphase
          phi = phi0 + 2d0*pi*dble(iphase-1)/Nphase
          call set_initial_condition(zpsi_store,zHO_store, &
                                     zpsi_CTEF, zHO_CTEF, phi, norm)

          call propagation(norm_CTEF_t,Ekin_CTEF_t)
!          if(myrank == 0)write(*,*)"norm",norm_CTEF_t(0),norm
          norm_CTEF_l = norm_CTEF_l + norm_CTEF_t*exp(-zI*phi)*norm*zweight

        end do
        
      end do

      norm_CTEF_l = norm_CTEF_l/dble(Ntraj*Nphase)
      call MPI_ALLREDUCE(norm_CTEF_l,norm_CTEF,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myrank == 0)then
        open(21,file="CTEF_norm.out")
        do it = 0,Nt+1
          write(21,"(999e26.16e3)")dt*it,norm_CTEF(it)
        end do
        close(21)
      end if

    end subroutine CTEF_dynamics
!-----------------------------------------------------------------------------------------
    subroutine init_forward_backward_trajectries(zpsi_out,zHO_out,zweight)
      implicit none
      complex(8),intent(out) :: zpsi_out(Lsite,2), zHO_out(Lsite,2),zweight
      integer :: i
      real(8) :: x1,x2,p1,p2

!! sub-system
!      do i = 1,Lsite
!        zpsi_out(i,:) = exp(zI*2d0*pi*(i-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite))
!      end do
      zpsi_out = 0d0
      zpsi_out(1,:) = 1d0

!check
!      p1 = 0d0
!      do i = 1,10000
!        call correlated_gaussian_random_number(x1,x2)
!        p1 = p1 + exp(-x1**2-x2**2)*exp(+0.5d0*x1**2+0.5d0*x2**2+0.5d0*(x1-x2)**2)*(2d0*pi/sqrt(3d0))/pi
!      end do
!      write(*,*)"p1=",p1/10000d0
!      stop

!! bath-system
      do i = 1,Lsite
        call correlated_gaussian_random_number(x1,x2)
        call correlated_gaussian_random_number(p1,p2)
        zHO_out(i,1) = x1 + zI * p1
        zHO_out(i,2) = x2 + zI * p2
      end do

      zweight = 1d0
      do i = 1, Lsite
        zweight = zweight * (2d0*pi/sqrt(3d0)/pi)**2*exp( &
          -0.5d0*abs(zHO_out(i,1))**2 -0.5d0*abs(zHO_out(i,2))**2 &
          +0.5d0*abs(zHO_out(i,1))**2 +0.5d0*abs(zHO_out(i,2))**2 &
          +0.5d0*abs(zHO_out(i,1)-zHO_out(i,2))**2 &
          )
      end do

!!! bath-system
!      do i = 1,Lsite
!        call gaussian_random_number(x1,p1)
!        call gaussian_random_number(x2,p2)
!        zHO_out(i,1) = x1 + zI * p1
!        zHO_out(i,2) = x2 + zI * p2
!      end do
!      zHO_out = zHO_out *sqrt(0.5d0)
!
!      zweight = 1d0
!      do i = 1, Lsite
!        zweight = zweight * (2d0/2d0)**2*exp( &
!          -0.5d0*abs(zHO_out(i,1))**2 -0.5d0*abs(zHO_out(i,2))**2 &
!          +1.0d0*abs(zHO_out(i,1))**2 +1.0d0*abs(zHO_out(i,2))**2 &
!          )
!      end do

    end subroutine init_forward_backward_trajectries
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
    subroutine propagation(norm_CTEF_t,Ekin_CTEF_t)
      implicit none
      real(8),intent(out) :: norm_CTEF_t(0:Nt+1),Ekin_CTEF_t(0:Nt+1)
      integer :: it

      call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(0))

      zHO_dot_CTEF = 0d0
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)

      do it = 0,Nt-1

        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
        call calc_norm(zpsi_CTEF,zHO_CTEF,norm_CTEF_t(it+1))
        
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
          zX_HO_CTEF(:,i,j) = sqrt(1d0/(2d0*mass*omega0)) &
            *(conjg(zHO_in(:,i)) + zHO_in(:,j)) &
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
      zEc_CTEF = -gamma*sqrt(2d0*mass*omega0) * zEc_CTEF 

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
      call dt_evolve_bath(zHO_inout,dt*0.5d0)
      zpsi_t = zpsi_inout
      zHO_t = zHO_inout


      call dt_evolve_elec(zpsi_inout,dt*0.5d0)
      call dt_evolve_bath(zHO_inout,dt*0.5d0)

      do iscf = 1, Nscf_pred_corr
        call refine_effective_hamiltonian(zpsi_inout,zHO_inout,zHO_dot_inout)
        zHO_inout = zHO_t
        zpsi_inout = zpsi_t

        call dt_evolve_elec(zpsi_inout,dt*0.5d0)
        call dt_evolve_bath(zHO_inout,dt*0.5d0)

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
    subroutine dt_evolve_bath(zHO_inout,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_inout(Lsite,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Nexp_Taylor = 6
      complex(8) :: zHO_t(Lsite,2),zhHO_t(Lsite,2)
      integer :: iexp
      complex(8) :: zfact

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_CTEF

      zHO_t = zHO_inout
      zfact = 1d0
      do iexp = 1,Nexp_Taylor
        zfact = zfact*(-zI*dt_t)/iexp
        zhHO_t(:,1) = zHb_eff_CTEF(1,1)*zHO_t(:,1) + zHb_eff_CTEF(1,2)*zHO_t(:,2)
        zhHO_t(:,2) = zHb_eff_CTEF(2,1)*zHO_t(:,1) + zHb_eff_CTEF(2,2)*zHO_t(:,2)

        zHO_inout = zHO_inout + zfact*zhHO_t
        zHO_t = zhHO_t
      end do

      zHO_inout = zHO_inout -zI*0.5d0*dt_t*zF_HO_CTEF

    end subroutine dt_evolve_bath
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
      real(8) :: c0

      c0 = -gamma*sqrt(2d0*mass*omega0)
      call hs_zpsi(zpsi_in,zhs_psi_t)

      zhpsi_t(:,1) = zSb_CTEF(1,1)*zhs_psi_t(:,1) + zSb_CTEF(1,2)*zhs_psi_t(:,2)
      zhpsi_t(:,2) = zSb_CTEF(2,1)*zhs_psi_t(:,1) + zSb_CTEF(2,2)*zhs_psi_t(:,2)

      zhpsi_t(:,1) = zhpsi_t(:,1) &
        +c0*zX_HO_CTEF(:,1,1)*zhs_psi_t(:,1) +c0*zX_HO_CTEF(:,1,2)*zhs_psi_t(:,2)
      zhpsi_t(:,2) = zhpsi_t(:,2) &
        +c0*zX_HO_CTEF(:,2,1)*zhs_psi_t(:,2) +c0*zX_HO_CTEF(:,2,2)*zhs_psi_t(:,2)

      zhpsi_t(:,1) = zhpsi_t(:,1) &
        -zDb_CTEF(1,1)*zpsi_in(:,1) -zDb_CTEF(1,2)*zpsi_in(:,2)
      zhpsi_t(:,2) = zhpsi_t(:,2) &
        -zDb_CTEF(2,1)*zpsi_in(:,1) -zDb_CTEF(2,2)*zpsi_in(:,2)

      zhpsi_out(:,1) = zSb_inv_CTEF(1,1)*zhpsi_t(:,1) + zSb_inv_CTEF(1,2)*zhpsi_t(:,2)
      zhpsi_out(:,2) = zSb_inv_CTEF(2,1)*zhpsi_t(:,1) + zSb_inv_CTEF(2,2)*zhpsi_t(:,2)

    end subroutine heff_zpsi
!-----------------------------------------------------------------------------------------
end module CTEF_mod
