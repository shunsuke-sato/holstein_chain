!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_dynamics
  use global_variables
  implicit none
  integer :: itraj,it
  integer :: isite
  real(8) :: pop(0:Nt+1,Lsite),pop_l(0:Nt+1,Lsite)
  real(8) :: Ekin_PBME(0:Nt+1),Ekin_PBME_l(0:Nt+1)
  real(8) :: Etot_PBME(0:Nt+1),Etot_PBME_l(0:Nt+1)
  complex(8) :: zpop0(Lsite), zEkin0_PBME, zEtot0_PBME
  integer :: icout  = 0
  real(8) :: ss

  pop_l = 0d0
  Ekin_PBME_l = 0d0
  Etot_PBME_l = 0d0
  do itraj = 1,Ntraj


    call PBME_initial_condition
    if(myrank == 0)then
      ss = sum(x_m**2+p_m**2)
      write(41,"(I9,999e16.6e3)")itraj,ss,2**(Lsite+1)*exp(-ss)*(ss-0.5d0*dble(Lsite)),abs(zweight0),real(zweight0)
    end if
    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0 .and. itraj > icout*(Ntraj/100))then
      write(*,"(I5,A)")icout,"% done"
      icout = icout + 1
    end if



    select case(PBME_flag)
    case('original','org','consistent')
      call PBE_population(zpop0)
      call PBE_Ekin(zEkin0_PBME,zEtot0_PBME)
    case('modified')
      call PBE_population(zpop0)
      call PBE_M_Ekin(zEkin0_PBME,zEtot0_PBME)
    case default
      stop 'Invalid PBME_flag'
    end select
    pop_l(0,:) = pop_l(0,:)  + zpop0(:)*zweight0
    Ekin_PBME_l(0) = Ekin_PBME_l(0)  + zEkin0_PBME*zweight0
    Etot_PBME_l(0) = Etot_PBME_l(0)  + zEtot0_PBME*zweight0


    do it = 0,Nt

      select case(PBME_flag)
      case('original','org')
        call PBME_dt_evolve_quantum !_traceless
      case('consistent')
        call PBME_C_dt_evolve_quantum !_traceless
      case('modified')
        call PBME_dt_evolve_quantum !_traceless
      case default
          stop 'Invalid PBME_flag'
      end select


      select case(PBME_flag)
      case('original','org','consistent')
        call PBE_population(zpop0)
        call PBE_Ekin(zEkin0_PBME,zEtot0_PBME)
      case('modified')
        call PBE_population(zpop0)
        call PBE_M_Ekin(zEkin0_PBME,zEtot0_PBME)
      case default
          stop 'Invalid PBME_flag'
      end select
      pop_l(it+1,:) = pop_l(it+1,:)  + zpop0(:)*zweight0
      Ekin_PBME_l(it+1) = Ekin_PBME_l(it+1)  + zEkin0_PBME*zweight0
      Etot_PBME_l(it+1) = Etot_PBME_l(it+1)  + zEtot0_PBME*zweight0


    end do

  end do

  call MPI_ALLREDUCE(pop_l,pop,(Nt+2)*Lsite,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  pop = pop/Ntraj
  call MPI_ALLREDUCE(Ekin_PBME_l,Ekin_PBME,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Ekin_PBME = Ekin_PBME/Ntraj
  call MPI_ALLREDUCE(Etot_PBME_l,Etot_PBME,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Etot_PBME = Etot_PBME/Ntraj

  if(myrank == 0)then
    open(21,file="PBME_population.out")
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,sum(pop(it,:)),pop(it,:)
    end do
    close(21)
  end if
  if(myrank == Nprocs-1)then
    open(21,file="PBME_Ekin.out")
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,Ekin_PBME(it),Etot_PBME(it)
    end do
    close(21)
  end if


end subroutine PBME_dynamics
!===============================================================================
subroutine PBE_population(zpop0)
  use global_variables
  implicit none
  complex(8) :: zpop0(Lsite)
  integer :: isite

  zpop0 = 0d0
  do isite = 1,Lsite
     zpop0(isite) = zpop0(isite) &
         + 0.5d0*(x_m(isite)**2 + p_m(isite)**2 -1d0) 
  end do

end subroutine PBE_population
!===============================================================================
subroutine PBE_M_population(zpop0)
  use global_variables
  implicit none
  complex(8) :: zpop0(Lsite)
  real(8) :: x2,weight
  integer :: isite

  x2 = sum(x_m**2 + p_m**2)
  weight = 2**(Lsite+1)*exp(-x2)

  x2 = 0d0
  zpop0 = 0d0
  do isite = 1,Lsite
     zpop0(isite) = zpop0(isite) &
         + (x_m(isite)**2 + p_m(isite)**2 -0.5d0) 
     x2 = x2 + zpop0(isite)
!     zpop0(isite) = zpop0(isite) &
!         + (x_m(isite)**2 + p_m(isite)**2 -0.5d0) * X_HO(isite)
  end do
  zpop0 = zpop0*weight
  x2 = x2 * weight

  do isite = 1,Lsite
     zpop0(isite) = x2*0.5d0*(X_HO(isite)**2 + V_HO(isite)**2 )
  end do


end subroutine PBE_M_population
!===============================================================================
subroutine PBE_Ekin(zEkin0,zEtot0)
  use global_variables
  implicit none
  complex(8) :: zEkin0,zEtot0
  complex(8) :: zdensity_matrix(Lsite,Lsite)
  complex(8) :: zdensity_matrix0(Lsite,Lsite)
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j
  complex(8) :: zs

  do i = 1,Lsite
    do j = 1,Lsite
      zdensity_matrix(i,j) = x_m(i)*x_m(j) + p_m(i)*p_m(j) &
        +zI*(x_m(i)*p_m(j) - x_m(j)*p_m(i) )
    end do
  end do

  do i = 1,Lsite
    zdensity_matrix(i,i) = zdensity_matrix(i,i) - 1d0
  end do
  zdensity_matrix0 = 0.5d0*zdensity_matrix


  zdensity_matrix = zdensity_matrix0*Hmat_kin
  zEkin0 = sum(zdensity_matrix)

!  return

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do


  zdensity_matrix = zdensity_matrix0 *Htot_t
  zEtot0 = sum(zdensity_matrix)

  zs = 0d0
  do i = 1,Lsite
     zs = zs + zdensity_matrix0(i,i)
  end do
  zEtot0 = zEtot0 + sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass) !*zs
  
!  Ekin0 = zweight0*sum(density_matrix*Hmat_kin)
!
!  Ekin0 = 0d0
!  do i = 1,Lsite
!    do j = 1,Lsite
!      Ekin0 = Ekin0 + zdensity_matrix(j,i)*zweight_m(j,i)
!    end do
!  end do

end subroutine PBE_Ekin
!===============================================================================
subroutine PBE_M_Ekin(zEkin0,zEtot0)
  use global_variables
  implicit none
  complex(8) :: zEkin0,zEtot0
  complex(8) :: zdensity_matrix(Lsite,Lsite)
  complex(8) :: zdensity_matrix0(Lsite,Lsite)
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j
  complex(8) :: zs
  real(8) :: x2,weight

  x2 = sum(x_m**2 + p_m**2)
  weight = 2**(Lsite+1)*exp(-x2)

  do i = 1,Lsite
    do j = 1,Lsite
      zdensity_matrix(i,j) = x_m(i)*x_m(j) + p_m(i)*p_m(j) !&
!        +zI*(x_m(i)*p_m(j) - x_m(j)*p_m(i) )
    end do
  end do

  do i = 1,Lsite
    zdensity_matrix(i,i) = zdensity_matrix(i,i) - 0.5d0
  end do
  zdensity_matrix0 = weight*zdensity_matrix


  zdensity_matrix = zdensity_matrix0*Hmat_kin
  zEkin0 = sum(zdensity_matrix)

!  return
  do i = 1,Lsite
    do j = 1,Lsite
      zdensity_matrix(i,j) = x_m(i)*x_m(j) + p_m(i)*p_m(j) &
        +zI*(x_m(i)*p_m(j) - x_m(j)*p_m(i) )
    end do
  end do

  do i = 1,Lsite
    zdensity_matrix(i,i) = zdensity_matrix(i,i) - 1d0
  end do
  zdensity_matrix0 = 0.5d0*zdensity_matrix


  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do


  zdensity_matrix = zdensity_matrix0 *Htot_t
  zEtot0 = sum(zdensity_matrix)

  zEtot0 = zEtot0 + sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass) !*zs
  Etot_mod = zEtot0
!  Ekin0 = zweight0*sum(density_matrix*Hmat_kin)
!
!  Ekin0 = 0d0
!  do i = 1,Lsite
!    do j = 1,Lsite
!      Ekin0 = Ekin0 + zdensity_matrix(j,i)*zweight_m(j,i)
!    end do
!  end do

end subroutine PBE_M_Ekin
!===============================================================================
subroutine PBE_M2_Ekin(zEkin0,zEtot0)
  use global_variables
  implicit none
  complex(8) :: zEkin0,zEtot0
  complex(8) :: zdensity_matrix(Lsite,Lsite)
  complex(8) :: zdensity_matrix0(Lsite,Lsite)
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j
  complex(8) :: zs
  real(8) :: x2,weight

  x2 = sum(x_m**2 + p_m**2)
  weight = 2**(Lsite+1)*exp(-x2)

  do i = 1,Lsite
    do j = 1,Lsite
      zdensity_matrix(i,j) = x_m(i)*x_m(j) + p_m(i)*p_m(j) !&
!        +zI*(x_m(i)*p_m(j) - x_m(j)*p_m(i) )
    end do
  end do

  do i = 1,Lsite
    zdensity_matrix(i,i) = zdensity_matrix(i,i) - 0.5d0
  end do
  zdensity_matrix0 = weight*zdensity_matrix


  zdensity_matrix = zdensity_matrix0*Hmat_kin
  zEkin0 = sum(zdensity_matrix)

!  return

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do


  zdensity_matrix = zdensity_matrix0 *Htot_t
  zEtot0 = sum(zdensity_matrix)

  zs = 0d0
  do i = 1,Lsite
     zs = zs + zdensity_matrix0(i,i)
  end do
  zEtot0 = zEtot0 + sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass) !*zs
  Etot_mod = zEtot0
!  Ekin0 = zweight0*sum(density_matrix*Hmat_kin)
!
!  Ekin0 = 0d0
!  do i = 1,Lsite
!    do j = 1,Lsite
!      Ekin0 = Ekin0 + zdensity_matrix(j,i)*zweight_m(j,i)
!    end do
!  end do

end subroutine PBE_M2_Ekin

