!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine allocation_of_global_arrays
  use global_variables
  implicit none
  integer :: i,j

  allocate(zC(Lsite),Hmat_kin(Lsite,Lsite),Hmat_coup(Lsite,Lsite),Hmat_tot(Lsite,Lsite))
  allocate(X_HO(Lsite),V_HO(Lsite),F_HO(Lsite),X_HO_ini(Lsite),V_HO_ini(Lsite))
  allocate(X_HO_new(Lsite),X_HO_old(Lsite))
  allocate(Ekin(0:Nt+1),Eph(0:Nt+1),Ecoup(0:Nt+1))
  allocate(Ekin_l(0:Nt+1),Eph_l(0:Nt+1),Ecoup_l(0:Nt+1))


  Hmat_kin = 0d0
  do i =1,Lsite
    j=i+1
    if(j>Lsite)j=j-Lsite
    Hmat_kin(i,j) = -t0
    j=i-1
    if(j<1)j=j+Lsite
    Hmat_kin(i,j) = -t0
  end do

end subroutine allocation_of_global_arrays
