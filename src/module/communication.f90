!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module communication
  use mpi
  use parallel
  implicit none

  private


  public :: comm_bcast

  interface comm_bcast
    module procedure comm_bcast_integer
    module procedure comm_bcast_integer_1d
    module procedure comm_bcast_integer_2d
    module procedure comm_bcast_integer_3d
    module procedure comm_bcast_real8
    module procedure comm_bcast_real8_1d
    module procedure comm_bcast_real8_2d
    module procedure comm_bcast_real8_3d
    module procedure comm_bcast_real8_4d
    module procedure comm_bcast_complex8
    module procedure comm_bcast_complex8_1d
    module procedure comm_bcast_complex8_2d
    module procedure comm_bcast_complex8_3d
    module procedure comm_bcast_complex8_4d
    module procedure comm_bcast_character
    module procedure comm_bcast_character_1d
    module procedure comm_bcast_logical
 end interface comm_bcast

contains
!-------------------------------------------------------------------------------
  function int_switch(if_true, int_true, int_false) result(int_result)
    implicit none
    logical,intent(in) :: if_true
    integer,intent(in) :: int_true, int_false
    integer :: int_result

    if(if_true)then
      int_result = int_true
    else
      int_result = int_false
   end if
  end function int_switch
!-------------------------------------------------------------------------------
  include "include/comm_bcast.f90"
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end module communication
