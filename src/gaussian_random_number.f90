!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine gaussian_random_number(x1,x2)
  use luxury
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  integer :: len = 4
  real :: rvec4(4)
  real(8) :: rvec(2),tmp
  integer(8) :: int1, int2

!  call random_number(r1)
!  call random_number(r2)
  CALL RANLUX (rvec4, len)
!  write(*,"(999e36.26e3)")rvec4(1)
!  write(*,"(999e36.26e3)")dble(rvec4(1))
!  write(*,"(999e36.26e3)")rvec4(2)
!  write(*,"(999e36.26e3)")dble(rvec4(2))
!  write(*,"(999e36.26e3)")rvec4(3)
!  write(*,"(999e36.26e3)")dble(rvec4(3))
!  write(*,"(999e36.26e3)")rvec4(4)
!  write(*,"(999e36.26e3)")dble(rvec4(4))
!
  int1 = aint(rvec4(1)*1d7)*10000000 + aint(rvec4(2)*1d7)
  int2 = aint(rvec4(3)*1d7)*10000000 + aint(rvec4(4)*1d7)
  rvec(1) =dble(int1)*1d-14
  rvec(2) =dble(int2)*1d-14

!  write(*,"(999e36.26e3)")rvec(1)
!  write(*,"(999e36.26e3)")rvec(2)

  stop
  if(rvec(1) == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(rvec(1)))
    x1 = tmp*cos(2d0*pi*rvec(2))
    x2 = tmp*sin(2d0*pi*rvec(2))
  end if
  return
end subroutine gaussian_random_number
