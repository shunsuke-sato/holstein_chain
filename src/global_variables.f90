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

! Control parameter
  character(8) :: calc_mode

! Parameters of model
  integer :: Lsite
  real(8) :: t0,gamma,mass,omega0,Tph
  complex(8),allocatable :: zC(:),zCp(:)
  real(8),allocatable :: Hmat_kin(:,:),Hmat_coup(:,:),Hmat_tot(:,:)
  real(8),allocatable :: X_HO(:),V_HO(:),F_HO(:),X_HO_new(:),X_HO_old(:)
  real(8),allocatable :: Xp_HO(:),Vp_HO(:),Fp_HO(:),Xp_HO_new(:),Xp_HO_old(:)

! Multi-trajectory
  integer :: Ntraj
  real(8),allocatable :: X_Ho_ini(:),V_HO_ini(:)

! Time propagation
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: Ekin(:),Eph(:),Ecoup(:),norm_t(:)
  real(8),allocatable :: Ekin_l(:),Eph_l(:),Ecoup_l(:),norm_t_l(:)

! Generalized Quantum master equation
  complex(8),allocatable :: zrho_DM(:,:,:)
  complex(8),allocatable :: zK_full(:,:,:,:,:),zK1(:,:,:,:,:),zK2(:,:,:,:,:),zK3(:,:,:,:,:)
  real(8) :: beta_KB

! PBME
  character(len=20) :: PBME_flag = 'original'
  real(8),allocatable :: x_m(:),p_m(:)
  complex(8),allocatable :: zweight_m(:,:)
  complex(8) :: zweight0
  real(8) :: x2_mean = -100d0,Etot_mod
!  real(8) :: x2_max
  

! FBTS
  character(len=20) :: FBTS_flag = 'original'
  real(8),allocatable :: x_n(:),p_n(:)


! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

! I/O
  integer,parameter :: nfile_full_kernel = 50
  character(50),parameter :: file_full_kernel = 'full_kernel.out'

! Timer 
  real(8) :: timer_start,timer_end


end module global_variables
