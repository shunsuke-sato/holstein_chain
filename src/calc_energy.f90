!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine calc_energy(Ekin_t,Eph_t,Ecoup_t)
  use global_variables
  implicit none
  real(8),intent(out) :: Ekin_t,Eph_t,Ecoup_t
  complex(8) :: zhC(Lsite)
  integer :: i

  zhC(:) = matmul(Hmat_kin,zC)
  Ekin_t = sum( conjg(zC)*zhC )
  Eph_t = sum(0.5d0*V_HO**2/mass + 0.5d0*X_HO**2*omega0**2*mass)
  Ecoup_t = -gamma*sqrt(2d0*mass*omega0)*sum(X_HO*abs(zC)**2)
  

end subroutine calc_energy
