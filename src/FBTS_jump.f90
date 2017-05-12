!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_jump(it)
  use global_variables
  implicit none
  integer,intent(in) :: it

  select case(FBTS_jump_flag)
  case('regular')
    call jump_regular
  case('exponential')
    call jump_exponential
  case default
    stop 'Invalid FBTS_flag'
  end select

  contains
!=========================================================================================
    subroutine jump_regular
      implicit none
      integer :: step,isite
      complex(8) :: z_m(Lsite),z_n(Lsite)
      complex(8) :: zs
      real(8) :: xx,pp

      step = aint(FBTS_jump_period/dt)
      if(mod(it+1,step) /= 0)return

      z_m = x_m + zI*p_m ; z_n = x_n + zI*p_n

      do isite = 1,Lsite
        call gaussian_random_number(xx,pp)
        x_m(isite) = sqrt(0.5d0)*xx; p_m(isite) = sqrt(0.5d0)*pp
        call gaussian_random_number(xx,pp)
        x_n(isite) = sqrt(0.5d0)*xx; p_n(isite) = sqrt(0.5d0)*pp
      end do

      zs = sum(conjg(z_m)*(x_m + zI*p_m)) * sum(z_n*conjg(x_n + zI*p_n))

      zweight0 = zweight0 * zs

    end subroutine jump_regular
!=========================================================================================
    subroutine jump_exponential
      implicit none
      integer :: step,isite
      complex(8) :: z_m(Lsite),z_n(Lsite)
      complex(8) :: zs
      real(8) :: xx,pp
      integer :: len = 1
      real(8) :: rvec(1)

      call ranlux_double(rvec,len)
      if(rvec(1) > dt/FBTS_jump_period)return

      z_m = x_m + zI*p_m ; z_n = x_n + zI*p_n

      do isite = 1,Lsite
        call gaussian_random_number(xx,pp)
        x_m(isite) = sqrt(0.5d0)*xx; p_m(isite) = sqrt(0.5d0)*pp
        call gaussian_random_number(xx,pp)
        x_n(isite) = sqrt(0.5d0)*xx; p_n(isite) = sqrt(0.5d0)*pp
      end do

      zs = sum(conjg(z_m)*(x_m + zI*p_m)) * sum(z_n*conjg(x_n + zI*p_n))

      zweight0 = zweight0 * zs

    end subroutine jump_exponential
!=========================================================================================
end subroutine FBTS_jump
