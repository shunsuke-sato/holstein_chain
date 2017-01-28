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
  allocate(X_HO(Lsite),V_HO(Lsite),X_HO_ini(Lsite,Ntraj),V_HO_ini(Lsite,Ntraj))
  allocate(X_HO_new(Lsite),V_HO_new(Lsite),X_HO_old(Lsite),V_HO_old(Lsite))
  allocate(F_HO(Lsite),F_HO_new(Lsite),F_HO_old(Lsite))

  Hmat_kin = 0d0
  do i =1,Lsite
    do j =i+1,Lsite
      if(abs(i-j) == 1)then
        Hmat_kin(i,j) = -t0
        Hmat_kin(j,i) = -t0
      end if
    end do
  end do

end subroutine allocation_of_global_arrays
