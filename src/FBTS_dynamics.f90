!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_dynamics
  use global_variables
  implicit none
  integer :: itraj,it
  integer :: isite
  real(8) :: pop(0:Nt+1,Lsite),pop_l(0:Nt+1,Lsite)
  real(8) :: Ekin_FBTS(0:Nt+1),Ekin_FBTS_l(0:Nt+1)
  real(8) :: Etot_FBTS(0:Nt+1),Etot_FBTS_l(0:Nt+1)
  complex(8) :: zpop0(Lsite),zEkin0_FBTS,zEtot0_FBTS
  integer :: icout  = 0
  real(8) :: Etot0

  pop_l = 0d0
  Ekin_FBTS_l = 0d0
  Etot_FBTS_l = 0d0

  do itraj = 1,Ntraj


    call FBTS_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0 .and. itraj > icout*(Ntraj/100))then
      write(*,"(I5,A)")icout,"% done"
      icout = icout + 1
    end if

    call FBTS_population(zpop0)
    pop_l(0,:) = pop_l(0,:)  + zpop0(:)*zweight0
    call FBTS_Ekin(zEkin0_FBTS,zEtot0_FBTS)
    Ekin_FBTS_l(0) = Ekin_FBTS_l(0)  + zEkin0_FBTS*zweight0
    Etot_FBTS_l(0) = Etot_FBTS_l(0)  + zEtot0_FBTS*zweight0

    if(itraj == 1)then
      open(30,file="traj01.out")
      call FBTS_calc_ene(Etot0)
      write(30,"(999e26.16e3)")0d0,Etot0
    end if

    do it = 0,Nt

      call FBTS_dt_evolve !_traceless
      call FBTS_population(zpop0)
      pop_l(it+1,:) = pop_l(it+1,:)  + zpop0(:)*zweight0
      call FBTS_Ekin(zEkin0_FBTS,zEtot0_FBTS)
      Ekin_FBTS_l(it+1) = Ekin_FBTS_l(it+1)  + zEkin0_FBTS*zweight0
      Etot_FBTS_l(it+1) = Etot_FBTS_l(it+1)  + zEtot0_FBTS*zweight0

    if(itraj == 1)then
      call FBTS_calc_ene(Etot0)
      write(30,"(999e26.16e3)")dt*(it+1),Etot0
    end if

    end do

    if(itraj == 1)then
      close(30)
    end if

  end do

  call MPI_ALLREDUCE(pop_l,pop,(Nt+2)*Lsite,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  pop = pop/Ntraj
  call MPI_ALLREDUCE(Ekin_FBTS_l,Ekin_FBTS,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Ekin_FBTS = Ekin_FBTS/Ntraj
  call MPI_ALLREDUCE(Etot_FBTS_l,Etot_FBTS,Nt+2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Etot_FBTS = Etot_FBTS/Ntraj

  if(myrank == 0)then
    open(21,file="FBTS_population.out")
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,sum(pop(it,:)),pop(it,:)
    end do
    close(21)
  end if
  if(myrank == Nprocs-1)then
    open(21,file="FBTS_Ekin.out")
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,Ekin_FBTS(it),Etot_FBTS(it)
    end do
    close(21)
  end if


end subroutine FBTS_dynamics
!===============================================================================
subroutine FBTS_population(zpop0)
  use global_variables
  implicit none
  complex(8) :: zpop0(Lsite)
  integer :: isite

  zpop0 = 0d0
  do isite = 1,Lsite
    zpop0(isite) = zpop0(isite) + (x_m(isite) - zI * p_m(isite)) &
      *(x_n(isite) + zI * p_n(isite)) 

  end do

end subroutine FBTS_population
!===============================================================================
subroutine FBTS_Ekin(zEkin0,zEtot0)
  use global_variables
  implicit none
  complex(8) :: zEkin0,zEtot0
  complex(8) :: zdensity_matrix(Lsite,Lsite)
  complex(8) :: zdensity_matrix0(Lsite,Lsite)
  real(8) :: n(Lsite),Htot_t(Lsite,Lsite)
  integer :: i,j
  real(8) :: ss
  complex(8) :: zs

  do i = 1,Lsite
    do j = 1,Lsite
      zdensity_matrix(i,j) = (x_m(i) - zI * p_m(i)) * (x_n(j) + zI * p_n(j))
    end do
  end do

  zdensity_matrix0 = zdensity_matrix


  zdensity_matrix = zdensity_matrix0*Hmat_kin
  zEkin0 = sum(zdensity_matrix)

!  return

  Htot_t = Hmat_kin
  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) - gamma*sqrt(2d0*mass*omega0)*X_HO(i)
  end do

  ss = sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)
  ss = ss - sum(  -gamma*sqrt(2d0*mass*omega0)*X_HO(:) ) ! traceless

  do i = 1,Lsite
    Htot_t(i,i) = Htot_t(i,i) + ss 
  end do


  zdensity_matrix = zdensity_matrix0 *Htot_t
  zEtot0 = sum(zdensity_matrix)
  return

!  zs = 0d0
!  do i = 1,Lsite
!    zs = zs + zdensity_matrix0(i,i)
!  end do
!  zEtot0 = zEtot0 + sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)*zs
!  zEtot0 = zEtot0 - sum(- gamma*sqrt(2d0*mass*omega0)*X_HO(:))*zs ! traceless


!  zEtot0 =  sum(0.5d0*mass*V_HO**2 + 0.5d0*X_HO**2*omega0**2*mass)*zs


!  Ekin0 = zweight0*sum(density_matrix*Hmat_kin)
!
!  Ekin0 = 0d0
!  do i = 1,Lsite
!    do j = 1,Lsite
!      Ekin0 = Ekin0 + zdensity_matrix(j,i)*zweight_m(j,i)
!    end do
!  end do

end subroutine FBTS_Ekin
