!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module global_variables
  implicit none
! Mathematical parameters
!  real(8),parameter :: pi=3.14159265358979323846d0
  real(8),parameter :: pi=4d0*atan(1d0)
  complex(8),parameter :: zI=(0d0,1d0)

! Parameters of model
  integer :: Lsite
  real(8) :: t0,gamma,mass,omega0,Tph
  complex(8),allocatable :: zC(:)
  real(8),allocatable :: Hmat_kin(:,:),Hmat_coup(:,:),Hmat_tot(:,:)
  real(8),allocatable :: X_HO(:),V_HO(:),X_HO_new(:),V_HO_new(:),X_HO_old(:),V_HO_old(:)
  real(8),allocatable :: F_HO(:),F_HO_new(:),F_HO_old(:)

! Multi-trajectory
  integer :: Ntraj
  real(8),allocatable :: X_Ho_ini(:,:),V_HO_ini(:,:)

! Time propagation
  integer :: Nt
  real(8) :: dt

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

end module global_variables
