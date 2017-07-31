!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module parallel
  use mpi
  implicit none

  private
! MPI global
  integer, public :: comm_group_global, &
                     comm_id_global, &
                     comm_nproc_global
  logical, public :: if_root_global
                     
! OMP
  integer, public :: nthread_omp

  public :: init_parallel, &
            fin_parallel,  &
            error_finalize

contains
!-------------------------------------------------------------------------------
  subroutine init_parallel
    implicit none
    integer :: ierr
!$ integer :: omp_get_max_threads  

    call MPI_init(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_nproc_global,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,comm_id_global,ierr)

    comm_group_global = MPI_COMM_WORLD

    if(comm_id_global == 0)then
       if_root_global = .true.
    else
       if_root_global = .false.
    end if

    nthread_omp = 1
!$  nthread_omp=omp_get_max_threads()

  end subroutine init_parallel
!-------------------------------------------------------------------------------
  subroutine fin_parallel
    implicit none
    integer :: ierr

    call MPI_Finalize(ierr)

  end subroutine fin_parallel
!-------------------------------------------------------------------------------
  subroutine error_finalize(message)
    implicit none
    character(*),intent(in) :: message
    integer :: ierr

    if(if_root_global)write(*,"(A)")message
    call MPI_Finalize(ierr)
    stop

  end subroutine error_finalize
!-------------------------------------------------------------------------------
end module parallel
