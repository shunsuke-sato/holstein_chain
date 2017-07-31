!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_model_parameters
  use global_variables
  use communication
  implicit none


  if(myrank == 0)then
    read(*,*)calc_mode
    read(*,*)Lsite
    read(*,*)t0
    read(*,*)omega0
    read(*,*)gamma
    read(*,*)mass
    read(*,*)Tph
    read(*,*)Ntraj
    read(*,*)Tprop
    read(*,*)dt
  end if

  call comm_bcast(calc_mode)
  call comm_bcast(Lsite)
  call comm_bcast(t0)
  call comm_bcast(omega0)
  call comm_bcast(gamma)
  call comm_bcast(mass)
  call comm_bcast(Tph)
  call comm_bcast(Ntraj)
  call comm_bcast(Tprop)
  call comm_bcast(dt)

  Nt = aint(Tprop/dt)+1


end subroutine set_model_parameters
!==============================================
subroutine mean_population
  use global_variables
  real(8) :: x2,y,ymin,ymax
  real(8) :: fact

  ymin = 1d0
  ymax = 2d0

  fact = 0.5d0*sqrt(exp(1d0)/4d0)**Lsite
  if(myrank == 0)write(*,*)"fact=",fact

  if(yimn*exp(-ymin) < fact)write(*,*)"something wrong..."
  do 
    y = ymax
    if(y*exp(-y) < fact )exit
    ymax = ymax*2d0
  end do


  do
    if(abs(ymax-ymin) < 1d-10)exit
    y = 0.5d0*(ymax + ymin)
    if(y*exp(-y) < fact)then
      ymax = y
    else
      ymin = y
    end if
  end do
  
  if(myrank == 0)write(*,*)"mean y= ",y
  x2 = y + dble(Lsite)/2d0
  if(myrank == 0)write(*,*)"mean x2= ",x2
  fact = 2**(Lsite+1)*exp(-x2)*(x2-dble(Lsite)/2d0)
  if(myrank == 0)write(*,*)"mean pop. = ",fact

  x2_mean = x2

end subroutine mean_population
