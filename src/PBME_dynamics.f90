!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_dynamics
  use global_variables
  implicit none
  integer :: itraj,it

  do itraj = 1,Ntraj


    call PBME_initial_conditions
    if(mod(itraj,Nprocs) /= myrank)cycle

    do it = 0,Nt

      call PBME_dt_evolve

    end do

  end do

end subroutine PBME_dynamics
