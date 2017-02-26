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
  integer :: icout  = 0
  

  pop_l = 0d0
  do itraj = 1,Ntraj


    call PBME_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0 .and. itraj > icout*(Ntraj/100))then
      write(*,"(I5,A)")icout,"% done"
      icout = icout + 1
    end if

! Start: calculate population
    it = -1
    do isite = 1,Lsite
      pop_l(it+1,isite) = pop_l(it+1,isite) &
        + zweight_m(1,1)*0.5d0*(x_m(isite)**2 + p_m(isite)**2 -1d0) 
    end do
! End: calculate population

    do it = 0,Nt

      call PBME_dt_evolve

! Start: calculate population
      do isite = 1,Lsite
        pop_l(it+1,isite) = pop_l(it+1,isite) &
          + zweight_m(1,1)*0.5d0*(x_m(isite)**2 + p_m(isite)**2 -1d0) 
      end do
! End: calculate population



    end do

  end do

  call MPI_ALLREDUCE(pop_l,pop,(Nt+2)*Lsite,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  pop = pop/Ntraj

  if(myrank == 0)then
    open(21,file="PBME_population.out")
    do it = 0,Nt+1
      write(21,"(999e26.16e3)")dt*it,pop(it,:)
    end do
    close(21)
  end if


end subroutine PBME_dynamics
