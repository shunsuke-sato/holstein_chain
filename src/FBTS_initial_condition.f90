!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_initial_condition
  use global_variables
  implicit none
  integer :: isite,jsite,i
  real(8) :: xx,pp
  real(8) :: n_m(Lsite),n_n(Lsite),n(Lsite)

  call set_initial_conditions_ph

  do isite = 1,Lsite
    call gaussian_random_number(xx,pp)
    x_m(isite) = sqrt(0.5d0)*xx; p_m(isite) = sqrt(0.5d0)*pp
    call gaussian_random_number(xx,pp)
    x_n(isite) = sqrt(0.5d0)*xx; p_n(isite) = sqrt(0.5d0)*pp
  end do

! Initial condition for phonon-bath
  do i = 1,Lsite
    n_m(i) = 0.5d0*(x_m(i)**2 + p_m(i)**2)
    n_n(i) = 0.5d0*(x_n(i)**2 + p_n(i)**2)
!    n(i) = n_m(i) + n_n(i) ! with trace
    n(i) = n_m(i) + n_n(i) - 1.0d0 ! traceless
    F_HO(i) = gamma*sqrt(2d0*mass*omega0)*n(i) - omega0**2*mass*X_HO(i) 
  end do
  X_HO_old = X_HO - V_HO*dt +0.5d0*F_HO/mass*dt**2


  
  do isite = 1,Lsite
    do jsite = 1,Lsite
      zweight_m(isite,jsite) = (x_m(isite) + zI * p_m(isite) )&
                              *(x_n(jsite) - zI * p_n(jsite) )
    end do
  end do

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

!  zweight0 = conjg(zweight0)
end subroutine FBTS_initial_condition
