!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine set_model_parameters
  use global_variables
  implicit none

  Lsite = 12
  t0 = 1d0
  omega0 = 1d0
  gamma = sqrt(0.4d0)
  mass = 1d0

  Tph = -1d0

  Ntraj = 1000000

  dt = 0.04d0 !0.01d0 !0.04d0
!  Nt = aint(30d0/dt)+1
!  Nt = aint(25d0/dt)+1
!  Nt = aint(5d0/dt)+1
  Nt = aint(25d0/dt)+1

!'MTEF', 'GQME_K'

!  calc_mode = 'MTEF'
  calc_mode = 'PTEF'
!  calc_mode = 'GQME_K'
!  calc_mode = 'GQME_T'
!  calc_mode = 'PBME'; PBME_flag = 'modified' ! 'original', 'consisten', 'modified'
!  x2_max = 5d0*dble(Lsite)
! calc_mode = 'FBTS'; FBTS_flag = 'modified' ! 'original', 'consisten', 'modified'
!  calc_mode = 'FBTS_K'

!  call mean_population

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
