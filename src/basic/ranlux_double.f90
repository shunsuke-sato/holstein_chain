!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine ranlux_double(rvec,len)
  use luxury
  implicit none
  integer,intent(in) :: len
  real(8),intent(out) :: rvec(len)
  integer :: len4 = 2
  real :: rvec4(2)
  integer(8) :: int1,i


  do i = 1,len

    CALL RANLUX (rvec4, len4)

    int1 = aint(rvec4(1)*1d6)*10000000 + aint(rvec4(2)*1d7)
    rvec(i) =dble(int1)*1d-13

  end do

end subroutine ranlux_double
