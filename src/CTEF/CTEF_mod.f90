!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module CTEF_mod
  use global_variables
  implicit none

  private

  integer,parameter :: Nphase = 2
  complex(8),allocatable :: zpsi_CTEF(:,:),zHO_CTEF(:,:)
  complex(8),allocatable :: zHO_dot_CTEF(:,:)

  complex(8) :: zSs_CTEF(2,2), zDs_CTEF(2,2)
  complex(8) :: zSb_CTEF(2,2), zDb_CTEF(2,2)
  complex(8) :: zEs_CTEF(2,2), zEc_CTEF(2,2), zEb_CTEF(2,2)
  complex(8),allocatable :: zX_HO_CTEF(:,:,:)

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
      allocate(zX_HO_CTEF(Lsite,2,2))

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
      integer :: itraj, iphase
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
      if(myrank == 0)write(*,"(A,2x,99e26.16e3)")"norm",norm_CTEF(0)

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

!        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF,zHO_dot_CTEF)
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
      complex(8) :: zs
      integer :: i,j

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


    end subroutine refine_effective_hamiltonian
!-----------------------------------------------------------------------------------------
    subroutine hs_zpsi(zpsi_in,zhpsi_out)
      implicit none
      integer :: i
      complex(8),intent(in) :: zpsi_in(Lsite,2)
      complex(8),intent(out) :: zhpsi_out(Lsite,2)


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
end module CTEF_mod
