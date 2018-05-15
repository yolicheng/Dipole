!***********************************************************************
!     ****   muller   ****
! Use Muller's method for rooting finding, three initial points are 
! used to contruct a parabola, and take the intersection of the x-axis with the parabola as the next approximation.

! Since Muller only deals with one parameter, we need to use do loop for all the components 

!  ****** Input Variables ******
!  n            dimension of the input point 
!  x1,x2, x3    of dimension n, the three initial points 
!  y(3)         the function to find the root, f(x(i)) = 0,i=1,...,n 
!  

! From Wikipedia 

!  ****** Output Variables ******


!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine muller( n, x1, x2, x3, y,  x)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)     ::  n   
real(kind=dp), intent(in)      :: x1(n), x2(n), x3(n), y(3)
real(kind=dp), intent(out)     :: x(n)   
 
real(kind=dp) ::  dnrm2 
! Local Variable
integer :: i
real(kind=dp)  :: xtemp(3), ytemp(3), xkm3, xkm2, xkm1, ykm3,ykm2,ykm1, & 
                  dy12, dy23, dy13, fk12, fk23, fk13,                  & 
                  fk123, dx, det, denomi, omega, xk, yk , xm
  
  ytemp = y 
  ! fk12 is the divided difference, fk12 = (y2-y1) / (x2-x1)
  
  print*, 'Input of muller,  y'
  print*, y
  print*, x1 
  print*, x2
  print*, x3 
  print*; read*
  
  do i = 1, n, 1
    xtemp = (/x1(i), x2(i), x3(i)/) 
    ytemp = y 
    
    xm = dnrm2(3, xtemp, 1) 
    print*, 'xtemp', xtemp 
    print*, 'ytemp', ytemp
    print*, 'xm',xm; print*; read*
    if(xm < 1.d-13) then 
      x(i) = 0.d0 
      cycle
    endif 
    
    ! the recurrent process
    do 
      
      xkm3 = xtemp(1); xkm2 = xtemp(2); xkm1 = xtemp(3)  
      ykm3 = ytemp(1); ykm2 = ytemp(2); ykm1 = ytemp(3)
   
      dy12 = ykm2 - ykm1
      dy23 = ykm3 - ykm2 
      dy13 = ykm3 - ykm1
    
      fk12 = dy12 / (xkm2 - xkm1)
      fk23 = dy23 / (xkm3 - xkm2)
      fk13 = dy13 / (xkm3 - xkm1)
    
      fk123 = (fk23-fk12) / (xkm3 - xkm1)
      omega = fk12 + fk13 + fk23
    
      print*, 'dy12, dy23, dy13', dy12, dy23, dy13
      print*, 'fk12, fk23, fk13, fk123', fk12, fk23, fk13, fk123
      print*, 'omega', omega
      print*; read*
      
      ! the old three points gives a parabola approximation 
      ! of the function 
    
!      xk = (/xkm3, xkm2, xkm1/)
    
      ! here, we need to check if the guess is in real number, by check det   > 0, 
      det = omega**2 - 4.d0*ykm1*fk123
      
      
      if(det < 0.d0) then 
        print*, 'Det < 0! obtain complex number! '
        print*, det
        print*; read* 
      endif  
    
      denomi = dsqrt(det) 
      if(omega < 0.d0)  then 
        denomi = omega - denomi
      else 
        denomi = omega + denomi
      endif     
    
      ! the new guess 
      xk = xkm1 - 2.d0*ykm1 / denomi 
      dx = xk - xkm1 
      
      print*, 'xk,denomi, dx ', xk, denomi, dx
      print*, 'omega, fk123 ', omega, fk123   
      print*; read*
      ! the function value at the new guess
      yk = ykm1  + omega*dx + fk123*dx*dx
      
      print*, 'check the new yk and old y(3)'
      print*, yk, ytemp; print*; read*
      
      if (dabs(yk) < 1.d-13)  exit   
      
      ! other wise, we update xtemp and ytemp, and continue
      xtemp = (/xkm2, xkm1, xk/)
      ytemp = (/ykm2, ykm1, yk/)
      
      print*, 'xtemp, ytemp '
      print*, xtemp; print*, ytemp
      cycle 
       
    enddo 
    
    ! get the new guess, update x(i)
    x(i) = xk   
    
  end do
 

  return  
end subroutine muller

