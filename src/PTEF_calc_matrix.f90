!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PTEF_calc_matrix(zSm,zHm,zHk,zHph,zHcoup,zDm)
  use global_variables
  implicit none
  complex(8),intent(out) :: zSm(2,2),zHm(2,2)
  complex(8),intent(out) :: zHk(2,2),zHph(2,2),zHcoup(2,2),zDm(2,2)
  complex(8) :: zhC(Lsite),zCt_t(Lsite)
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)
  complex(8) :: zdt_HO(Lsite),zpdt_HO(Lsite)
  complex(8) :: zexponent, zMm(2,2), zovl_s,zovl_b
  integer :: i

  z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
  zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
  zdt_HO = (V_HO + zI*F_HO/(omega0*mass))*sqrt(mass*omega0/2d0)
  zpdt_HO = (Vp_HO + zI*Fp_HO/(omega0*mass))*sqrt(mass*omega0/2d0)

  zexponent = -0.5d0*sum(abs(z_HO - zp_HO)**2) -zI*aimag(sum(z_HO*conjg(zp_HO)))

  zSm(1,1) = 1d0 !sum(abs(zC)**2) !1d0
  zSm(2,2) = 1d0 !sum(abs(zCp)**2) !1d0
  zovl_s = sum(conjg(zC)*zCp)
  zovl_b =  exp(zexponent)
  zSm(1,2) = zovl_s * zovl_b
  zSm(2,1) = conjg(zSm(1,2))
!  write(*,*)zSm(1,2),zSm(2,1)

! Hk
  zhC(:) = matmul(Hmat_kin,zC)
  zHk(1,1) = sum(conjg(zC)*zhC)
  zhC(:) = matmul(Hmat_kin,zCp)
  zHk(2,2) = sum(conjg(zCp)*zhC)
  zHk(1,2) = sum(conjg(zC)*zhC)*zovl_b
  zHk(2,1) = conjg (zHk(1,2) )

! Hph
  zHph(1,1) = omega0*sum(abs(z_HO)**2) + dble(Lsite)/2d0
  zHph(2,2) = omega0*sum(abs(zp_HO)**2) + dble(Lsite)/2d0
  zHph(1,2) = omega0*sum(conjg(z_HO)*zp_HO) + dble(Lsite)/2d0
  zHph(1,2) = zHph(1,2)*zovl_s*zovl_b
  zHph(2,1) = conjg(zHph(1,2))


! Hcoup
  zHcoup(1,1) = -gamma*sum((conjg(z_HO) + z_HO)*abs(zC)**2)
  zHcoup(2,2) = -gamma*sum( (conjg(zp_HO) + zp_HO)*abs(zCp)**2 )
  zHcoup(1,2) = -gamma*sum( (conjg(z_HO) + zp_HO )*conjg(zC)*zCp )*zovl_b
  zHcoup(2,1) = conjg(zHcoup(1,2))
!  zHcoup =  0d0 ! test


! Dm
  zCt_t(:) = zC(:)
  call hpsi(zCt_t,zhC)
  zDm(1,1) = real(sum(conjg(zC)*zhC)) + aimag(sum(conjg(zdt_HO)*z_HO))
!  zDm(1,1) = real(sum(conjg(zC)*zhC)) + zI*2d0*zI*aimag(sum(conjg(z_HO)*zdt_HO))

  zDm(2,1) = sum(conjg(zCp)*zhC)*conjg(zovl_b) &
    +zi*conjg(zovl_s*zovl_b)*sum(conjg(zp_HO)*zdt_HO - real(zdt_HO*conjg(z_HO)))
!  zDm(2,1) = sum(conjg(zCp)*zhC)*conjg(zovl_b) &
!    +zi*conjg(zovl_s*zovl_b)*sum(zdt_HO*conjg(zp_HO) - conjg(zdt_HO)*z_HO  )


  zCt_t(:) = zCp(:)
  call hpsi_pair(zCt_t,zhC)
  zDm(2,2) = real(sum(conjg(zCp)*zhC)) + aimag(sum(conjg(zpdt_HO)*zp_HO))
!  zDm(1,1) = real(sum(conjg(zC)*zhC)) + zI*2d0 * zI*aimag(sum(conjg(zp_HO)*zpdt_HO))

  zDm(1,2) = sum(conjg(zC)*zhC)*zovl_b &
    +zi*zovl_s*zovl_b*sum(conjg(z_HO)*zpdt_HO - real(zpdt_HO*conjg(zp_HO)))
!  zDm(1,2) = sum(conjg(zC)*zhC)*zovl_b &
!    +zi*zovl_s*zovl_b*sum(zpdt_HO*conjg(z_HO) - conjg(zpdt_HO)*zp_HO   )


  zHm = zHk + zHph + zHcoup

! Hcoup


end subroutine PTEF_calc_matrix


!=====================================================================================
subroutine PTEF_calc_matrix_old(zSm,zVm,zVm_i_v22,zHm,zHm_v,zHk_v,zHph_v,zHcoup_v,zSm_v)
  use global_variables
  implicit none
  complex(8),intent(out) :: zSm(2,2),zVm(2,2),zVm_i_v22(2,2)
  complex(8),intent(out) :: zHm(2,2),zHm_v(2,2)
  complex(8),intent(out) :: zHk_v(2,2),zHph_v(2,2),zHcoup_v(2,2),zSm_v(2,2)
  complex(8) :: zhC(Lsite)
  complex(8) :: z_HO(Lsite),zp_HO(Lsite)
  complex(8) :: zexponent, zMm(2,2), zovl_s,zovl_b
  integer :: i

  z_HO = (X_HO + zI*V_HO/omega0)*sqrt(mass*omega0/2d0)
  zp_HO = (Xp_HO + zI*Vp_HO/omega0)*sqrt(mass*omega0/2d0)
  zexponent = -0.5d0*sum(abs(z_HO - zp_HO)**2) -zI*aimag(sum(z_HO*conjg(zp_HO)))

  zSm(1,1) = 1d0; zSm(2,2) = 1d0
  zovl_s = sum(conjg(zC)*zCp)
  zovl_b = exp(zexponent)
  zSm(1,2) = zovl_s * zovl_b 
  zSm(2,1) = conjg(zSm(1,2))


  zVm(1,1) = 1d0; zVm(1,2) = zSm(1,2)
  zVm(2,1) = 0d0; zVm(2,2) = sqrt(abs(1d0 - abs(zSm(1,2))**2))
  zVm_i_v22(1,1) = zVm(2,2); zVm_i_v22(1,2) = - zVm(1,2)
  zVm_i_v22(2,1) = 0d0     ; zVm_i_v22(2,2) =   zVm(1,1)

! Hk
  zhC(:) = matmul(Hmat_kin,zC)
  zHk_v(1,1) = sum(conjg(zC)*zhC)
  zhC(:) = matmul(Hmat_kin,zCp)
  zHk_v(2,2) = sum(conjg(zCp)*zhC)
  zHk_v(1,2) = sum(conjg(zC)*zhC)*zovl_b
  zHk_v(2,1) = conjg (zHk_v(1,2) )

! Hph
  zHph_v(1,1) = omega0*sum(abs(z_HO)**2) + dble(Lsite)/2d0
  zHph_v(2,2) = omega0*sum(abs(z_HO)**2) + dble(Lsite)/2d0
  zHph_v(1,2) = omega0*sum(conjg(z_HO)*zp_HO) + dble(Lsite)/2d0
  zHph_v(1,2) = zHph_v(1,2)*zovl_s*zovl_b
  zHph_v(2,1) = conjg(zHph_v(1,2))

! Hcoup
  zHcoup_v(1,1) = -gamma*sum(conjg(z_HO) + z_HO)
  zHcoup_v(2,2) = -gamma*sum(conjg(zp_HO) + zp_HO)
  zHcoup_v(1,2) = -gamma*sum(conjg(z_HO) + zp_HO)*zovl_s*zovl_b
  zHcoup_v(2,1) = conjg(zHcoup_v(1,2))

  zHm = zHk_v + zHph_v + zHcoup_v

  zMm = matmul(zHk_v,zVm_i_v22); zHk_v = matmul( transpose(conjg(zVm_i_v22)),zMm )
  zMm = matmul(zHph_v,zVm_i_v22); zHph_v = matmul( transpose(conjg(zVm_i_v22)),zMm )
  zMm = matmul(zHcoup_v,zVm_i_v22); zHcoup_v = matmul( transpose(conjg(zVm_i_v22)),zMm )
  zMm = matmul(zSm,zVm_i_v22); zSm_v = matmul( transpose(conjg(zVm_i_v22)),zMm )

  zHm_v = zHk_v + zHph_v + zHcoup_v
  zHm_v = zHm_v/abs(zVm(2,2))**2



end subroutine PTEF_calc_matrix_old
