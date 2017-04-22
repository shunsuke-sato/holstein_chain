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

subroutine calc_energy_pair(zEkin_t,zEph_t,zEcoup_t)
  use global_variables
  implicit none
  complex(8),intent(out) :: zEkin_t,zEph_t,zEcoup_t
  complex(8) :: zhCp(Lsite)
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)
  integer :: i

  z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
  zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)

  zhCp(:) = matmul(Hmat_kin,zCp)

  zEkin_t = sum( conjg(zC)*zhCp )
!  zEph_t = omega0*sum(conjg(z_HO)*zp_HO) + dble(Lsite)/2d0
!  zEph_t = zEph_t*sum(conjg(zC)*zCp)
  zEph_t = 1d0 !1d0*sum(conjg(zC)*zCp)
  zEcoup_t = -gamma*sum((conjg(z_HO) + zp_HO)*conjg(zC)*zCp)
  

end subroutine calc_energy_pair
