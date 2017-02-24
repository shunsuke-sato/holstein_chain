!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_thermal_ph_dist
  use global_variables
  implicit none
  integer :: i,itraj
  real(8) :: beta,fact,ss
  real(8) :: x1,x2

  if(Tph <= 0d0)then
    fact = 2d0/omega0
  else
    beta = 1d0/Tph
    fact = 2d0*tanh(beta*omega0/2d0)/omega0
  end if

  beta_KB = fact

  do i = 1,Lsite

     call gaussian_random_number(x1,x2)
     ss = x1/sqrt(fact) !Quantum distribution
     X_HO_ini(i)=ss/(omega0*sqrt(mass))
     ss = x2/sqrt(fact) !Quantum distribution
     V_HO_ini(i)=ss/sqrt(mass)

  end do

end subroutine set_thermal_ph_dist
