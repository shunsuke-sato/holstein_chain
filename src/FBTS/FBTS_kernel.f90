!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_kernel
  use global_variables
  implicit none
  integer :: itraj,it
  integer :: isite
  integer :: icout  = 0

  complex(8),allocatable :: zK1_l(:,:,:,:,:),zK3_l(:,:,:,:,:)
  complex(8) :: zfact(1,Lsite), zDM(Lsite,Lsite)
  integer :: i,j
  integer :: a1,a2,b1,b2

  call allocate_GQME_kernel

  allocate(zK1_l(Lsite,Lsite,1,Lsite,0:Nt),zK3_l(Lsite,Lsite,1,Lsite,0:Nt))



  zK1_l = 0d0; zK3_l = 0d0    
  do itraj = 1,Ntraj


    call FBTS_initial_condition
    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0 .and. itraj > icout*(Ntraj/100))then
      write(*,"(I5,A)")icout,"% done"
      icout = icout + 1
    end if


    do j = 1,Lsite
      zfact(1,j) = (x_m(j) + zI * p_m(j) )&
                  *(x_n(1) - zI * p_n(1) )
!      zfact(1,j) = zfact(1,j) &
!        * ( (X_HO(1) - X_HO(j))  +zI*0.5d0*beta_KB*(V_HO(1)+V_HO(j)) )

      zfact(1,j) = zfact(1,j) &
        * ( (X_HO(1) - X_HO(j))  -zI*0.5d0*beta_KB*(V_HO(1)+V_HO(j)) )

    end do

    do i = 1,Lsite
      do j = 1,Lsite
        zDM(i,j) = (x_m(i) - zI * p_m(i)) * (x_n(j) + zI * p_n(j))
      end do
    end do


    do b2 = 1,Lsite
      do a1=1,Lsite
        do a2=1,Lsite
          zK1_l(a1,a2,1,b2,0) = zK1_l(a1,a2,1,b2,0) &
            +gamma**2*2d0*mass*omega0*zDM(a2,a1) &
            *(X_HO(a1) - X_HO(a2))*zfact(1,b2)
          zK3_l(a1,a2,1,b2,0) = zK3_l(a1,a2,1,b2,0) &
            -gamma*sqrt(2d0*mass*omega0)*zDM(a2,a1) &
            *zfact(1,b2)
        end do
      end do
    end do

    do it = 0,Nt-1

      call FBTS_dt_evolve



      do i = 1,Lsite
         do j = 1,Lsite
            zDM(i,j) = (x_m(i) - zI * p_m(i)) * (x_n(j) + zI * p_n(j))
         end do
      end do

      do b2 = 1,Lsite
        do a1=1,Lsite
          do a2=1,Lsite
            zK1_l(a1,a2,1,b2,it + 1) = zK1_l(a1,a2,1,b2,it + 1) &
              +gamma**2*2d0*mass*omega0*zDM(a2,a1) &
              *(X_HO(a1) - X_HO(a2))*zfact(1,b2)
            zK3_l(a1,a2,1,b2,it + 1) = zK3_l(a1,a2,1,b2, it + 1) &
              -gamma*sqrt(2d0*mass*omega0)*zDM(a2,a1) &
              *zfact(1,b2)
          end do
        end do
      end do

    end do

  end do

  call MPI_ALLREDUCE(zK1_l,zK1,(Nt+1)*Lsite**3, &
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(zK3_l,zK3,(Nt+1)*Lsite**3, &
    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  zK1 = zK1/Ntraj
  zK3 = zK3/Ntraj

  call refine_GQME_kernel_K1K3  ! Inversion symmetry, Unitarity, ... etc.
  call evaluate_GQME_kernel_full

  if(myrank == 0)then
    open(nfile_full_kernel,file=trim(file_full_kernel),form='unformatted')
    write(nfile_full_kernel)zK_full,zK1,zK3
    close(nfile_full_kernel)
  end if


end subroutine FBTS_kernel
