!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_initial_condition
  use global_variables
  implicit none
  integer :: isite,jsite
  real(8) :: xx,pp

  call set_initial_conditions_ph

  do isite = 1,Lsite
    call gaussian_random_number(xx,pp)
    x_m(isite) = sqrt(0.5d0)*xx; p_m(isite) = sqrt(0.5d0)*pp
  end do
  
  do isite = 1,Lsite
    do jsite = 1,Lsite
      zweight_m(isite,jsite) = x_m(isite)*x_m(jsite) &
         + p_m(isite)*p_m(jsite) &
         + zI*(x_m(isite)*p_m(jsite) - p_m(isite)*x_m(jsite) )
      if(isite == jsite)zweight_m(isite,jsite) = zweight_m(isite,jsite) -0.5d0
    end do
  end do
  zweight_m = 2d0 * zweight_m


end subroutine PBME_initial_condition
