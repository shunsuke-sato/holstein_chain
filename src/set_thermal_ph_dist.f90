!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_thermal_ph_dist
  use global_variables
  implicit none
  integer :: i,itraj
  real(8) :: beta,fact,ss
  real(8) :: x1,x2

  if(Tph <= 0d0)then
    fact = 2d0/omega0
  else
    beta = 1d0/Tph
    fact = 2d0*tanh(beta*omega0/2d0)/omega0
  end if


  if(myrank == 0)then
    do itraj=1,Ntraj
      do i = 1,Lsite

        call normal_random_number(x1,x2)
        call random_number(x2)
        ss = x1/sqrt(fact) !Quantum distribution
        X_HO_ini(i,itraj)=ss*sin(2d0*pi*x2)/(omega0*sqrt(mass))
        V_HO_ini(i,itraj)=ss*cos(2d0*pi*x2)*sqrt(mass)

      end do
    end do
  end if

  call MPI_BCAST(X_HO_ini,Lsite*Ntraj,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(V_HO_ini,Lsite*Ntraj,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

end subroutine set_thermal_ph_dist


subroutine normal_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  real(8) :: r1,r2,tmp

  call random_number(r1)
  call random_number(r2)

  if(r1 == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(r1))
    x1 = tmp*cos(2d0*pi*r2)
    x2 = tmp*sin(2d0*pi*r2)
  end if
  return
end subroutine normal_random_number
