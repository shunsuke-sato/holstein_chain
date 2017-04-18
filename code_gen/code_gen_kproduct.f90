!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
program main
  implicit none
  integer,parameter :: Lsite = 12
  character(256) :: filename = 'kernel_product12.f90'
  character(256) :: routinename = 'kernel_product12'
  character(256) :: car(Lsite**2),ctmp,cnum(Lsite)
  integer :: ic,i
  integer :: mod_table(-Lsite:Lsite)
  integer :: a1,a2,b1,b2,c1,c2,a1t,a2t,c1t,c2t

  do i = -Lsite,Lsite
    mod_table(i) = mod(i+2*Lsite,Lsite)
  end do
  do i = 1,Lsite
    write(cnum(i),"(i0)")i
  end do

  car(1)='subroutine '//trim(routinename)//'(zK1,zK2,zK3)'
  car(2)='implicit none'
  write(ctmp,"(I2)")Lsite
  car(3)='complex(8),intent(in) :: zK1('//trim(ctmp)//','//trim(ctmp)//',1,'//trim(ctmp)//')'
  car(4)='complex(8),intent(in) :: zK2('//trim(ctmp)//','//trim(ctmp)//',1,'//trim(ctmp)//')'
  car(5)='complex(8),intent(in) :: zK3('//trim(ctmp)//','//trim(ctmp)//',1,'//trim(ctmp)//')'

  open(21,file=filename)
  do ic = 1,5
    write(21,"(A)")trim(car(ic))

  end do


  b1=1
  do b2=1,Lsite
    do a1 = 1,Lsite
      do a2 = 1,Lsite

        car(1)='zK3('//trim(cnum(a1))//','//trim(cnum(a2))//',1,'//trim(cnum(b2))//') = &'
        write(21,"(A)")trim(car(1))


        ic = 0
        do c2 = 1,Lsite
          do c1 = 1,Lsite
            ic = ic +1
            c1t=1
            c2t = mod_table(c2-c1) + 1
            a1t = mod_table(a1-c1) + 1
            a2t = mod_table(a2-c1) + 1


            car(ic)=&
              '+zKk1('//trim(cnum(a1t))//','//trim(cnum(a2t))//',1,'//trim(cnum(c2t))//')*'
            car(ic)= &
             trim(car(ic))//'zK2('//trim(cnum(c1))//','//trim(cnum(c2))//',1,'//trim(cnum(b2))//')'
            if(ic /= Lsite**2)car(ic) = trim(car(ic))//'&'

          end do
        end do

        do ic = 1,Lsite**2
          write(21,"(A)")trim(car(ic))
        end do

      end do
    end do
  end do



  car(1) = 'end subroutine '//trim(routinename)
  write(21,"(A)")trim(car(1))
  close(21)
end program main
