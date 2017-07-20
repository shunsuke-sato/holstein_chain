!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine correlated_gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  integer :: len = 3
  real(8) :: rvec(3),r1,r2,r3
  integer(8) :: int1, int2
  real(8) :: tmp

  do 
    call ranlux_double (rvec, len)
!    call random_number(r1)
!    call random_number(r2)
!    call random_number(r3)
!    rvec(1) = r1
!    rvec(2) = r2
!    rvec(3) = r3

    if(rvec(1) == 0d0)then
      x1 = 0d0
      x2 = 0d0
    else 
      tmp = sqrt(-2d0*log(rvec(1)))
      x1 = tmp*cos(2d0*pi*rvec(2))
      x2 = tmp*sin(2d0*pi*rvec(2))
    end if

    tmp = x1 -x2
    r1 = exp(-0.5d0*tmp**2)
    r2 = rvec(3)

    if(r2 < r1)exit
  end do


  return
end subroutine correlated_gaussian_random_number
