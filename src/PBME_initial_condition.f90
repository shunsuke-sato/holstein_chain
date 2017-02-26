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
         - zI*(x_m(isite)*p_m(jsite) - p_m(isite)*x_m(jsite) )
      if(isite == jsite)zweight_m(isite,jsite) = zweight_m(isite,jsite) -0.5d0
    end do
  end do
  zweight_m = 2d0 * zweight_m

!  zweight0 = zweight_m(1,1)
!  return

  zweight0 = 0d0
  do isite = 1,Lsite
    do jsite = 1,Lsite

      
!      zweight_m(isite,jsite) = zweight_m(isite,jsite)&
      zweight0 = zweight0 + zweight_m(isite,jsite)&
        *exp(zI*2d0*pi*(isite-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite)) &
        *exp(-zI*2d0*pi*(jsite-1)*dble(Lsite/2)/dble(Lsite))/sqrt(dble(Lsite))

    end do
  end do

end subroutine PBME_initial_condition
